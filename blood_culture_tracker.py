"""
Blood Culture Contamination Tracker core calculations.

The module provides transparent rule-based heuristics for blood-culture
adjudication, differential time-to-positivity (DTTP), contamination
surveillance, Wilson confidence intervals, and optional economic modeling.

Important: outputs are decision-support signals, not clinical diagnoses or
CDC/NHSN surveillance determinations.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple


class AdjudicationVerdict(str, Enum):
    TRUE_PATHOGEN = "Likely true bloodstream infection signal"
    PROBABLE_CONTAMINATION = "Probable blood-culture contamination signal"
    # Kept for backward API compatibility. DTTP supports CRBSI assessment; it
    # is not by itself a CDC/NHSN CLABSI surveillance criterion.
    CLABSI_SUSPECTED = "Catheter-related bloodstream infection (CRBSI) supported by DTTP"
    INDETERMINATE = "Indeterminate; clinical correlation and repeat culture may be needed"


TRUE_PATHOGENS: Set[str] = {
    "staphylococcus_aureus",
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "pseudomonas_aeruginosa",
    "streptococcus_pneumoniae",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "candida_albicans",
    "candida_glabrata",
    "streptococcus_pyogenes",
    "streptococcus_agalactiae",
    "neisseria_meningitidis",
    "bacteroides_fragilis",
    "listeria_monocytogenes",
    "salmonella_enterica",
    "acinetobacter_baumannii",
    "serratia_marcescens",
    "proteus_mirabilis",
}

COMMON_CONTAMINANTS: Set[str] = {
    "coagulase_negative_staphylococcus",
    "staphylococcus_epidermidis",
    "staphylococcus_hominis",
    "staphylococcus_capitis",
    "staphylococcus_warneri",
    "corynebacterium_species",
    "cutibacterium_acnes",
    "propionibacterium_acnes",
    "micrococcus_luteus",
    "bacillus_species",
    "viridans_group_streptococci",
    "aerococcus_viridans",
}

# Illustrative defaults retained for backward compatibility. These are not
# universal cost estimates; callers should replace them with local values.
DEFAULT_COST_PER_CONTAMINATION_USD: Dict[str, float] = {
    "excess_length_of_stay": 4500.0,
    "unnecessary_antimicrobials": 1200.0,
    "additional_lab_and_microbiology": 550.0,
    "echocardiogram_and_consults": 850.0,
}


def normalize_organism(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")
    aliases = {
        "coagulase_negative_staphylococci": "coagulase_negative_staphylococcus",
        "cons": "coagulase_negative_staphylococcus",
        "corynebacterium_spp": "corynebacterium_species",
        "bacillus_spp": "bacillus_species",
        "viridans_streptococci": "viridans_group_streptococci",
        "propionibacterium_acnes": "cutibacterium_acnes",
    }
    return aliases.get(text, text) or None


@dataclass
class CultureBottle:
    bottle_id: str
    site_type: str
    ttp_hours: Optional[float] = None
    is_positive: bool = False
    organism: Optional[str] = None


@dataclass
class BloodCultureSet:
    set_id: str
    patient_id: str
    collection_unit: str
    bottles: List[CultureBottle] = field(default_factory=list)
    clinical_fever_temp_c: Optional[float] = None
    has_indwelling_catheter: bool = False

    @property
    def total_bottles(self) -> int:
        return len(self.bottles)

    @property
    def positive_bottles(self) -> int:
        return sum(1 for bottle in self.bottles if bottle.is_positive)

    @property
    def primary_organism(self) -> Optional[str]:
        for bottle in self.bottles:
            if bottle.is_positive and bottle.organism:
                return normalize_organism(bottle.organism)
        return None

    @property
    def shortest_ttp_hours(self) -> Optional[float]:
        ttps = [
            bottle.ttp_hours
            for bottle in self.bottles
            if bottle.is_positive and bottle.ttp_hours is not None
        ]
        return min(ttps) if ttps else None


class CultureAdjudicationEngine:
    """Rule-based decision-support classifier for a blood-culture set."""

    @classmethod
    def adjudicate_set(cls, culture_set: BloodCultureSet) -> Dict[str, Any]:
        reasons: List[str] = []
        score = 0

        org = culture_set.primary_organism
        if not org or culture_set.positive_bottles == 0:
            return {
                "set_id": culture_set.set_id,
                "patient_id": culture_set.patient_id,
                "organism": org,
                "verdict": "Negative culture / no positive bottles represented",
                "is_contamination": False,
                "confidence_score": 0,
                "shortest_ttp_hours": None,
                "dttp_hours": None,
                "bottles_summary": f"{culture_set.positive_bottles}/{culture_set.total_bottles} positive",
                "reasoning": ["No positive culture bottle with an organism was represented."],
                "clinical_recommendation": "Interpret with the full clinical and microbiology context.",
            }

        is_known_pathogen = org in TRUE_PATHOGENS
        is_known_contaminant = org in COMMON_CONTAMINANTS

        if is_known_pathogen:
            score += 3
            reasons.append(f"Organism '{org}' is in the configured pathogen list (+3).")
        elif is_known_contaminant:
            score -= 2
            reasons.append(f"Organism '{org}' is in the configured common-commensal list (-2).")
        else:
            reasons.append(f"Organism '{org}' is not classified by the built-in organism lists.")

        min_ttp = culture_set.shortest_ttp_hours
        if min_ttp is not None:
            if min_ttp < 14.0:
                score += 2
                reasons.append(f"Rapid time to positivity ({min_ttp:.1f} h < 14 h) increases the true-BSI signal (+2).")
            elif min_ttp > 36.0 and is_known_contaminant:
                score -= 2
                reasons.append(f"Delayed time to positivity ({min_ttp:.1f} h > 36 h) with a common commensal decreases the true-BSI signal (-2).")

        pos_ratio = culture_set.positive_bottles / max(culture_set.total_bottles, 1)
        if pos_ratio >= 0.75 and culture_set.total_bottles >= 2:
            score += 2
            reasons.append(
                f"High bottle concordance ({culture_set.positive_bottles}/{culture_set.total_bottles} positive) (+2)."
            )
        elif pos_ratio <= 0.25 and culture_set.total_bottles >= 4:
            score -= 2
            reasons.append(
                f"Low bottle concordance ({culture_set.positive_bottles}/{culture_set.total_bottles} positive) (-2)."
            )

        # DTTP compares matching-organism cultures. Positive DTTP means the
        # central-line bottle became positive earlier than the peripheral one.
        matching = [
            bottle
            for bottle in culture_set.bottles
            if bottle.is_positive
            and bottle.ttp_hours is not None
            and normalize_organism(bottle.organism) == org
        ]
        central_ttps = [
            bottle.ttp_hours for bottle in matching if bottle.site_type == "central_line"
        ]
        peripheral_ttps = [
            bottle.ttp_hours for bottle in matching if bottle.site_type == "peripheral"
        ]

        dttp_hours: Optional[float] = None
        catheter_source_supported = False
        if central_ttps and peripheral_ttps:
            dttp_hours = round(min(peripheral_ttps) - min(central_ttps), 2)
            if dttp_hours >= 2.0:
                catheter_source_supported = True
                score += 3
                reasons.append(
                    f"Central-line culture became positive {dttp_hours:.2f} h earlier than the paired peripheral culture; DTTP >= 2 h supports a catheter source (+3)."
                )
            elif dttp_hours <= -2.0:
                reasons.append(
                    f"Peripheral culture became positive {-dttp_hours:.2f} h earlier than the central-line culture; this does not support a catheter source by DTTP."
                )
            else:
                reasons.append(
                    f"Paired DTTP was {dttp_hours:.2f} h, below the 2 h catheter-source threshold."
                )

        if catheter_source_supported:
            verdict = AdjudicationVerdict.CLABSI_SUSPECTED
            is_contamination = False
            action = (
                "DTTP supports a catheter source. Correlate clinically and apply local CRBSI diagnostic guidance; "
                "do not treat this output as an NHSN CLABSI determination."
            )
        elif score >= 2:
            verdict = AdjudicationVerdict.TRUE_PATHOGEN
            is_contamination = False
            action = "The rule set favors true bloodstream infection; correlate with clinical findings and laboratory context."
        elif score <= -2:
            verdict = AdjudicationVerdict.PROBABLE_CONTAMINATION
            is_contamination = True
            action = "The rule set favors contamination; review repeat cultures, collection details, and the clinical context before changing therapy."
        else:
            verdict = AdjudicationVerdict.INDETERMINATE
            is_contamination = False
            action = "The rule set is indeterminate; use clinical assessment, repeat cultures when appropriate, and local policy."

        return {
            "set_id": culture_set.set_id,
            "patient_id": culture_set.patient_id,
            "organism": org,
            "verdict": verdict.value,
            "is_contamination": is_contamination,
            "confidence_score": score,
            "shortest_ttp_hours": min_ttp,
            "dttp_hours": dttp_hours,
            "bottles_summary": f"{culture_set.positive_bottles}/{culture_set.total_bottles} positive",
            "reasoning": reasons,
            "clinical_recommendation": action,
        }


class ContaminationSurveillanceEngine:
    """Contamination rates, Wilson intervals, and unit-level funnel limits."""

    DEFAULT_TARGET_PCT: float = 3.0
    # Backward-compatible name used by existing callers.
    CLSI_TARGET_PCT: float = DEFAULT_TARGET_PCT

    @classmethod
    def calculate_wilson_ci(
        cls,
        contaminated_count: int,
        total_sets: int,
        z_score: float = 1.96,
    ) -> Tuple[float, float]:
        if total_sets < 0 or contaminated_count < 0 or contaminated_count > total_sets:
            raise ValueError("Counts must satisfy 0 <= contaminated_count <= total_sets.")
        if total_sets == 0:
            return 0.0, 0.0
        if z_score <= 0:
            raise ValueError("z_score must be positive.")

        p = contaminated_count / total_sets
        z2 = z_score * z_score
        denom = 1.0 + z2 / total_sets
        center = (p + z2 / (2.0 * total_sets)) / denom
        margin = (
            z_score
            * math.sqrt((p * (1.0 - p) + z2 / (4.0 * total_sets)) / total_sets)
            / denom
        )
        return (
            round(max(0.0, (center - margin) * 100.0), 2),
            round(min(100.0, (center + margin) * 100.0), 2),
        )

    @classmethod
    def analyze_surveillance_data(
        cls,
        culture_records: List[Dict[str, Any]],
        target_pct: float = DEFAULT_TARGET_PCT,
    ) -> Dict[str, Any]:
        if not 0 <= target_pct <= 100:
            raise ValueError("target_pct must be between 0 and 100.")

        total_sets = len(culture_records)
        contaminated_sets = sum(
            1 for record in culture_records if bool(record.get("is_contamination", False))
        )
        rate_pct = contaminated_sets / total_sets * 100.0 if total_sets else 0.0
        ci_lower, ci_upper = cls.calculate_wilson_ci(contaminated_sets, total_sets)

        unit_buckets: Dict[str, Dict[str, int]] = {}
        for record in culture_records:
            unit = str(record.get("collection_unit") or "Unknown")
            aggregate = unit_buckets.setdefault(unit, {"total": 0, "contaminated": 0})
            aggregate["total"] += 1
            if bool(record.get("is_contamination", False)):
                aggregate["contaminated"] += 1

        overall_p = contaminated_sets / total_sets if total_sets else 0.0
        funnel_analysis: Dict[str, Dict[str, Any]] = {}
        for unit_name, counts in unit_buckets.items():
            n_u = counts["total"]
            k_u = counts["contaminated"]
            u_rate = k_u / n_u * 100.0 if n_u else 0.0
            se = math.sqrt(max(overall_p * (1.0 - overall_p), 1e-6) / max(n_u, 1)) * 100.0
            ucl_95 = min(100.0, rate_pct + 1.96 * se)
            ucl_997 = min(100.0, rate_pct + 3.0 * se)
            funnel_analysis[unit_name] = {
                "total_sets": n_u,
                "contaminated_sets": k_u,
                "rate_percentage": round(u_rate, 2),
                "upper_control_limit_95": round(ucl_95, 2),
                "upper_control_limit_997": round(ucl_997, 2),
                "outlier_95": u_rate > ucl_95,
                "outlier_997": u_rate > ucl_997,
            }

        meets_target = rate_pct <= target_pct
        return {
            "total_sets_evaluated": total_sets,
            "contaminated_sets_count": contaminated_sets,
            "overall_contamination_rate_pct": round(rate_pct, 2),
            "wilson_95ci_pct": [ci_lower, ci_upper],
            "configured_target_threshold_pct": target_pct,
            "clsi_target_threshold_pct": target_pct,
            "meets_configured_target": meets_target,
            "meets_clsi_standard": meets_target,
            "status_summary": "Within configured target" if meets_target else "Above configured target",
            "unit_funnel_analytics": funnel_analysis,
        }


class EconomicImpactEngine:
    """Optional cost model using caller-supplied or illustrative assumptions."""

    @classmethod
    def calculate_cost(
        cls,
        contaminated_count: int,
        custom_cost_weights: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        if contaminated_count < 0:
            raise ValueError("contaminated_count must be non-negative.")

        costs = dict(DEFAULT_COST_PER_CONTAMINATION_USD)
        if custom_cost_weights:
            for key, value in custom_cost_weights.items():
                value = float(value)
                if value < 0:
                    raise ValueError("Cost assumptions must be non-negative.")
                costs[key] = value

        per_case_total = sum(costs.values())
        total_cost = per_case_total * contaminated_count
        return {
            "contaminated_sets_count": contaminated_count,
            "cost_breakdown_per_contamination_usd": costs,
            "total_cost_per_case_usd": round(per_case_total, 2),
            "estimated_annual_excess_cost_usd": round(total_cost, 2),
            "potential_savings_50pct_reduction_usd": round(total_cost * 0.50, 2),
            "assumption_note": "Illustrative defaults; replace with local cost estimates before operational use.",
        }
