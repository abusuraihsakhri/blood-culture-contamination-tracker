#!/usr/bin/env python3
"""Backward-compatible DTTP and contamination-rate helpers.

This module delegates core calculations to blood_culture_tracker.py so the
repository has one source of truth for adjudication and statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from blood_culture_tracker import (
    BloodCultureSet,
    ContaminationSurveillanceEngine,
    CultureAdjudicationEngine,
    CultureBottle,
    EconomicImpactEngine,
    normalize_organism,
)


@dataclass
class CultureSet:
    set_id: str
    organism: Optional[str]
    ttp_hours: Optional[float]
    bottles_drawn: int
    bottles_positive: int
    central_only_positive: bool = False
    peripheral_ttp_hours: Optional[float] = None
    central_ttp_hours: Optional[float] = None


def adjudicate(cs: CultureSet) -> Dict[str, Any]:
    if cs.bottles_drawn <= 0 or not 0 <= cs.bottles_positive <= cs.bottles_drawn:
        raise ValueError("Bottle counts must satisfy 0 <= positive <= drawn and drawn > 0.")

    explicit = int(cs.peripheral_ttp_hours is not None) + int(cs.central_ttp_hours is not None)
    if explicit > cs.bottles_positive:
        raise ValueError("Paired TTP inputs exceed bottles_positive.")

    organism = normalize_organism(cs.organism)
    bottles: List[CultureBottle] = []

    if cs.peripheral_ttp_hours is not None:
        bottles.append(
            CultureBottle(
                f"{cs.set_id}-P",
                "peripheral",
                cs.peripheral_ttp_hours,
                True,
                organism,
            )
        )
    if cs.central_ttp_hours is not None:
        bottles.append(
            CultureBottle(
                f"{cs.set_id}-C",
                "central_line",
                cs.central_ttp_hours,
                True,
                organism,
            )
        )

    default_site = "central_line" if cs.central_only_positive else "peripheral"
    for index in range(cs.bottles_positive - len(bottles)):
        bottles.append(
            CultureBottle(
                f"{cs.set_id}-POS-{index + 1}",
                default_site,
                cs.ttp_hours,
                True,
                organism,
            )
        )
    for index in range(cs.bottles_drawn - len(bottles)):
        bottles.append(CultureBottle(f"{cs.set_id}-NEG-{index + 1}", "peripheral", None, False, None))

    result = CultureAdjudicationEngine.adjudicate_set(
        BloodCultureSet(cs.set_id, "PT-UNKNOWN", "Unknown", bottles)
    )

    if result["is_contamination"]:
        verdict = "probable_contamination"
    elif result["confidence_score"] >= 2:
        verdict = "likely_true_bacteremia"
    else:
        verdict = "indeterminate"

    return {
        "set_id": cs.set_id,
        "organism": organism,
        "central_peripheral_ttp_differential_h": result.get("dttp_hours"),
        "evidence_points": result["confidence_score"],
        "reasoning": result["reasoning"],
        "verdict": verdict,
        "detail": result["verdict"],
    }


def contamination_rate_stats(
    sets: List[Dict[str, Any]], target_pct: float = 3.0
) -> Dict[str, Any]:
    mapped = [
        {
            "is_contamination": bool(row.get("adjudicated_as_contamination")),
            "collection_unit": row.get("collection_unit", "UNKNOWN"),
        }
        for row in sets
    ]
    result = ContaminationSurveillanceEngine.analyze_surveillance_data(
        mapped, target_pct=target_pct
    )
    funnel = {
        unit: {
            "rate_pct": data["rate_percentage"],
            "n": data["total_sets"],
            "above_95pct_limit": data["outlier_95"],
            "limit_95pct": [0.0, data["upper_control_limit_95"]],
        }
        for unit, data in result["unit_funnel_analytics"].items()
    }
    return {
        "sets_analyzed": result["total_sets_evaluated"],
        "contaminated": result["contaminated_sets_count"],
        "overall_rate_pct": result["overall_contamination_rate_pct"],
        "wilson_95ci_pct": result["wilson_95ci_pct"],
        "target_pct": target_pct,
        "meets_configured_target": result["meets_configured_target"],
        "meets_clsi_target": result["meets_configured_target"],
        "by_collection_unit_funnel": funnel,
    }


def _wilson_ci_pct(k: int, n: int, z: float = 1.96) -> List[float]:
    lower, upper = ContaminationSurveillanceEngine.calculate_wilson_ci(k, n, z)
    return [lower, upper]


def economic_impact(
    contaminated_sets: int, costs: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    result = EconomicImpactEngine.calculate_cost(contaminated_sets, costs)
    return {
        "tool": "blood-culture-contamination-tracker",
        "contaminated_sets": contaminated_sets,
        "cost_breakdown_per_case_usd": result["cost_breakdown_per_contamination_usd"],
        "cost_per_contaminated_set_usd": result["total_cost_per_case_usd"],
        "estimated_total_cost_usd": result["estimated_annual_excess_cost_usd"],
        "avoidable_if_rate_halved_usd": result["potential_savings_50pct_reduction_usd"],
        "assumption_note": result["assumption_note"],
    }
