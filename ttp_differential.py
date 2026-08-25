#!/usr/bin/env python3
"""
Blood Culture Contamination Tracker — TTP Differential & Rate Analytics
Time-to-positivity differential algorithms (true bacteremia vs contaminant),
contamination-rate control limits (Wilson CI + funnel limits), and economic
impact quantification.

Zero-dependency. Author: Dr. Abu Suraih Sakhri. License: MIT.
"""
import argparse
import json
import math
import sys
from dataclasses import dataclass
from typing import Dict, Any, List, Optional


# Known skin-flora contaminants vs true pathogens
CONTAMINANT_ORGANISMS = {
    "coagulase_negative_staphylococcus", "corynebacterium", "cutibacterium_acnes",
    "micrococcus", "bacillus_species", "viridans_group_streptococci",
}
TRUE_PATHOGEN_EXAMPLES = {
    "staphylococcus_aureus", "escherichia_coli", "klebsiella_pneumoniae",
    "pseudomonas_aeruginosa", "streptococcus_pneumoniae", "enterococcus_faecalis",
}

# Published decision thresholds (hours)
TTP_TRUE_THRESHOLD_H = 15.0      # positivity < 15h: strong true bacteremia signal
TTP_CONTAMINANT_H = 40.0         # > 40h for skin flora: likely contaminant
COST_PER_CONTAMINATED = {        # published ranges (US$, 2010s-2020s literature)
    "extra_length_of_stay": 4800.0,
    "additional_antimicrobials": 1100.0,
    "additional_lab_testing": 600.0,
}


@dataclass
class CultureSet:
    set_id: str
    organism: Optional[str]
    ttp_hours: Optional[float]          # time to positivity of first bottle
    bottles_drawn: int
    bottles_positive: int
    central_only_positive: bool = False
    peripheral_ttp_hours: Optional[float] = None
    central_ttp_hours: Optional[float] = None


def adjudicate(cs: CultureSet) -> Dict[str, Any]:
    """Rule-based true-pathogen-vs-contaminant classification with reasoning."""
    reasons, points = [], 0
    is_contaminant_species = cs.organism in CONTAMINANT_ORGANISMS if cs.organism else False
    is_pathogen_species = cs.organism in TRUE_PATHOGEN_EXAMPLES if cs.organism else False

    if cs.ttp_hours is not None:
        if cs.ttp_hours < TTP_TRUE_THRESHOLD_H:
            points += 2; reasons.append(f"TTP {cs.ttp_hours:.1f}h < {TTP_TRUE_THRESHOLD_H:.0f}h (+2)")
        elif cs.ttp_hours > TTP_CONTAMINANT_H and is_contaminant_species:
            points -= 2; reasons.append(f"skin flora at TTP {cs.ttp_hours:.1f}h (>40h) (-2)")

    if is_pathogen_species:
        points += 2; reasons.append("obligate/typical pathogen species (+2)")
    if is_contaminant_species:
        points -= 1; reasons.append("recognized skin flora species (-1)")

    frac = cs.bottles_positive / max(cs.bottles_drawn, 1)
    if frac >= 0.5 and cs.bottles_drawn >= 2:
        points += 1; reasons.append(f">=50% bottles positive ({cs.bottles_positive}/{cs.bottles_drawn}) (+1)")
    elif frac <= 0.25 and cs.bottles_drawn >= 4:
        points -= 1; reasons.append("<25% bottles positive in full set (-1)")

    # bottle-pair differential logic
    differential = None
    if cs.peripheral_ttp_hours is not None and cs.central_ttp_hours is not None:
        differential = round(cs.central_ttp_hours - cs.peripheral_ttp_hours, 2)
        if abs(differential) >= 2.0:
            points += 2; reasons.append(
                f"|central-peripheral| differential {differential}h >= 2h suggests catheter source (+2)")
        elif abs(differential) < 1.0 and is_contaminant_species:
            points -= 1; reasons.append("concordant rapid positivity both sites; CoNS less concerning (-1)")
    elif cs.central_only_positive and is_contaminant_species:
        points += 1; reasons.append("CoNS from central bottles only: CLABSI possible (+1)")

    verdict = ("likely_true_bacteremia" if points >= 2 else
               "indeterminate" if points == 1 or points == -1 or points == 0 else
               "probable_contamination")
    return {
        "set_id": cs.set_id,
        "organism": cs.organism,
        "central_peripheral_ttp_differential_h": differential,
        "evidence_points": points,
        "reasoning": reasons,
        "verdict": verdict,
    }


def contamination_rate_stats(sets: List[Dict[str, Any]], target_pct: float = 3.0) -> Dict[str, Any]:
    """Overall + per-unit contamination rates with Wilson CIs and funnel limits."""
    n = len(sets)
    contaminated = sum(1 for s in sets if s.get("adjudicated_as_contamination"))
    rate = 100 * contaminated / n if n else 0.0
    ci = _wilson_ci_pct(contaminated, n)

    by_unit: Dict[str, Dict[str, int]] = {}
    for s in sets:
        unit = s.get("collection_unit", "UNKNOWN")
        agg = by_unit.setdefault(unit, {"total": 0, "contaminated": 0})
        agg["total"] += 1
        if s.get("adjudicated_as_contamination"):
            agg["contaminated"] += 1

    pbar = contaminated / n if n else 0.0
    funnel = {}
    for unit, agg in by_unit.items():
        nu = agg["total"]
        ku = agg["contaminated"]
        ru = 100 * ku / nu if nu else 0
        half = 1.96 * math.sqrt(max(pbar * (1 - pbar), 1e-9) / nu) * 100 if nu else 0
        funnel[unit] = {
            "rate_pct": round(ru, 2),
            "n": nu,
            "above_95pct_limit": ru > pbar * 100 + half,
            "limit_95pct": [round(max(pbar * 100 - half, 0.0), 2), round(pbar * 100 + half, 2)],
        }

    return {
        "sets_analyzed": n,
        "contaminated": contaminated,
        "overall_rate_pct": round(rate, 2),
        "wilson_95ci_pct": [round(x, 2) for x in ci],
        "target_pct": target_pct,
        "meets_clsi_target": rate < target_pct,
        "by_collection_unit_funnel": funnel,
    }


def _wilson_ci_pct(k: int, n: int, z: float = 1.96) -> List[float]:
    if not n:
        return [0.0, 0.0]
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return [max(0, (center - half)) * 100, min(1, center + half) * 100]


def economic_impact(contaminated_sets: int,
                    costs: Dict[str, float] = None) -> Dict[str, Any]:
    c = dict(COST_PER_CONTAMINATED)
    if costs:
        c.update(costs)
    per_case = sum(c.values())
    total = per_case * contaminated_sets
    return {
        "tool": "blood-culture-contamination-tracker/enriched",
        "contaminated_sets": contaminated_sets,
        "cost_breakdown_per_case_usd": c,
        "cost_per_contaminated_set_usd": round(per_case, 2),
        "estimated_total_cost_usd": round(total, 2),
        "avoidable_if_rate_halved_usd": round(total / 2, 2),
    }


if __name__ == "__main__":
    cases = [
        CultureSet("BC001", "staphylococcus_aureus", ttp_hours=9.5, bottles_drawn=2,
                   bottles_positive=2, central_ttp_hours=9.5, peripheral_ttp_hours=14.0),
        CultureSet("BC002", "coagulase_negative_staphylococcus", ttp_hours=44.0,
                   bottles_drawn=4, bottles_positive=1,
                   peripheral_ttp_hours=44.0, central_ttp_hours=None),
        CultureSet("BC003", "coagulase_negative_staphylococcus", ttp_hours=11.0,
                   bottles_drawn=2, bottles_positive=1, central_only_positive=True,
                   central_ttp_hours=11.0),
        CultureSet("BC004", "escherichia_coli", ttp_hours=7.8, bottles_drawn=2,
                   bottles_positive=2),
    ]
    for c in cases:
        print(json.dumps(adjudicate(c), indent=2))

    ledger = [
        {"id": "BC001", "adjudicated_as_contamination": False, "collection_unit": "ED"},
        {"id": "BC002", "adjudicated_as_contamination": True, "collection_unit": "ED"},
        {"id": "BC003", "adjudicated_as_contamination": False, "collection_unit": "ICU"},
        {"id": "BC004", "adjudicated_as_contamination": False, "collection_unit": "ICU"},
        {"id": "BC005", "adjudicated_as_contamination": True, "collection_unit": "WARD3"},
        {"id": "BC006", "adjudicated_as_contamination": True, "collection_unit": "WARD3"},
        {"id": "BC007", "adjudicated_as_contamination": False, "collection_unit": "WARD3"},
        {"id": "BC008", "adjudicated_as_contamination": False, "collection_unit": "ICU"},
    ]
    print(json.dumps(contamination_rate_stats(ledger), indent=2))
    print(json.dumps(economic_impact(3), indent=2))
