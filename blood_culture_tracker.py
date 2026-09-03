"""
Blood Culture Contamination Tracker & TTP Differential Engine
=============================================================
Clinical microbiology surveillance and diagnostic decision support adhering to
CLSI M47-A2 and CDC NHSN bloodstream infection standards:
- Time-To-Positivity (TTP) and Differential TTP (DTTP) for CLABSI diagnosis
- Organism pathogenicity matrix & bottle concordance adjudication
- Wilson score confidence intervals and Statistical Process Control (SPC) funnel limits
- Hospital economic impact & avoidable excess cost modeling
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple


class AdjudicationVerdict(str, Enum):
    TRUE_PATHOGEN = "Likely True Bacteremia / Sepsis"
    PROBABLE_CONTAMINATION = "Probable Skin Contaminant"
    CLABSI_SUSPECTED = "Suspected Catheter-Related Bloodstream Infection (CLABSI/CRBSI)"
    INDETERMINATE = "Indeterminate / Repeat Culture Recommended"


# CLSI / NHSN recognized obligate or frequent pathogens
TRUE_PATHOGENS: Set[str] = {
    "staphylococcus_aureus", "escherichia_coli", "klebsiella_pneumoniae",
    "pseudomonas_aeruginosa", "streptococcus_pneumoniae", "enterococcus_faecalis",
    "enterococcus_faecium", "candida_albicans", "candida_glabrata",
    "streptococcus_pyogenes", "streptococcus_agalactiae", "neisseria_meningitidis",
    "bacteroides_fragilis", "listeria_monocytogenes", "salmonella_enterica",
    "acinetobacter_baumannii", "serratia_marcescens", "proteus_mirabilis",
}

# Common skin commensals / typical contaminants
COMMON_CONTAMINANTS: Set[str] = {
    "coagulase_negative_staphylococcus", "staphylococcus_epidermidis",
    "staphylococcus_hominis", "staphylococcus_capitis", "staphylococcus_warneri",
    "corynebacterium_species", "cutibacterium_acnes", "propionibacterium_acnes",
    "micrococcus_luteus", "bacillus_species", "viridans_group_streptococci",
    "aerococcus_viridans",
}

# Economic costs per contaminated blood culture event (published clinical literature benchmark)
DEFAULT_COST_PER_CONTAMINATION_USD: Dict[str, float] = {
    "excess_length_of_stay": 4500.0,
    "unnecessary_antimicrobials": 1200.0,
    "additional_lab_and_microbiology": 550.0,
    "echocardiogram_and_consults": 850.0,
}


@dataclass
class CultureBottle:
    bottle_id: str
    site_type: str  # "peripheral" or "central_line"
    ttp_hours: Optional[float] = None
    is_positive: bool = False
    organism: Optional[str] = None


@dataclass
class BloodCultureSet:
    set_id: str
    patient_id: str
    collection_unit: str  # e.g., "ED", "ICU", "Ward_3B", "Phlebotomy"
    bottles: List[CultureBottle] = field(default_factory=list)
    clinical_fever_temp_c: Optional[float] = None
    has_indwelling_catheter: bool = False

    @property
    def total_bottles(self) -> int:
        return len(self.bottles)

    @property
    def positive_bottles(self) -> int:
        return sum(1 for b in self.bottles if b.is_positive)

    @property
    def primary_organism(self) -> Optional[str]:
        for b in self.bottles:
            if b.is_positive and b.organism:
                return b.organism.strip().lower().replace(" ", "_")
        return None

    @property
    def shortest_ttp_hours(self) -> Optional[float]:
        ttps = [b.ttp_hours for b in self.bottles if b.is_positive and b.ttp_hours is not None]
        return min(ttps) if ttps else None


# ==============================================================================
# 1. ADJUDICATION & DIFFERENTIAL TTP (DTTP) ENGINE
# ==============================================================================

class CultureAdjudicationEngine:
    """Classifies blood culture positivity as true bacteremia, contaminant, or CLABSI."""

    @classmethod
    def adjudicate_set(cls, culture_set: BloodCultureSet) -> Dict[str, Any]:
        """
        Adjudicate culture based on:
        1. Organism species taxonomy (Pathogen vs Commensal)
        2. Time to positivity (TTP) kinetics
        3. Differential TTP (DTTP) between central line and peripheral bottles
        4. Positive bottle concordance ratio
        """
        reasons: List[str] = []
        score = 0  # Positive = true bacteremia signal, Negative = contamination signal

        org = culture_set.primary_organism
        if not org or culture_set.positive_bottles == 0:
            return {
                "set_id": culture_set.set_id,
                "verdict": "Negative Culture (No Growth)",
                "confidence_score": 0,
                "is_contamination": False,
                "reasoning": ["All bottles negative for bacterial growth"],
                "clsi_action": "No antimicrobial escalation indicated",
            }

        is_known_pathogen = org in TRUE_PATHOGENS or any(p in org for p in TRUE_PATHOGENS)
        is_known_contaminant = org in COMMON_CONTAMINANTS or any(c in org for c in COMMON_CONTAMINANTS)

        if is_known_pathogen:
            score += 3
            reasons.append(f"Organism '{org}' is a clinically significant pathogen (+3)")
        elif is_known_contaminant:
            score -= 2
            reasons.append(f"Organism '{org}' is recognized common skin flora (-2)")

        # TTP Analysis
        min_ttp = culture_set.shortest_ttp_hours
        if min_ttp is not None:
            if min_ttp < 14.0:
                score += 2
                reasons.append(f"Rapid time-to-positivity ({min_ttp:.1f}h < 14h) suggests high bacterial load (+2)")
            elif min_ttp > 36.0 and is_known_contaminant:
                score -= 2
                reasons.append(f"Delayed time-to-positivity ({min_ttp:.1f}h > 36h) with skin commensal (-2)")

        # Bottle Concordance
        pos_ratio = culture_set.positive_bottles / max(culture_set.total_bottles, 1)
        if pos_ratio >= 0.75 and culture_set.total_bottles >= 2:
            score += 2
            reasons.append(f"High bottle concordance ({culture_set.positive_bottles}/{culture_set.total_bottles} bottles positive) (+2)")
        elif pos_ratio <= 0.25 and culture_set.total_bottles >= 4:
            score -= 2
            reasons.append(f"Low bottle concordance ({culture_set.positive_bottles}/{culture_set.total_bottles} bottles positive) (-2)")

        # Differential TTP for Catheter-Related Bloodstream Infection (DTTP)
        central_ttps = [b.ttp_hours for b in culture_set.bottles if b.is_positive and b.site_type == "central_line" and b.ttp_hours is not None]
        peripheral_ttps = [b.ttp_hours for b in culture_set.bottles if b.is_positive and b.site_type == "peripheral" and b.ttp_hours is not None]

        dttp_hours = None
        is_clabsi_suspect = False

        if central_ttps and peripheral_ttps:
            min_central = min(central_ttps)
            min_peripheral = min(peripheral_ttps)
            # DTTP = Peripheral TTP - Central TTP
            dttp_hours = round(min_peripheral - min_central, 2)
            if dttp_hours >= 2.0:
                is_clabsi_suspect = True
                score += 3
                reasons.append(f"Central line bottle turned positive {dttp_hours}h earlier than peripheral (DTTP >= 2.0h -> CLABSI criterion met) (+3)")
            elif dttp_hours <= -2.0:
                reasons.append(f"Peripheral bottle turned positive {-dttp_hours}h earlier than central (favors non-catheter source)")

        # Determine final verdict
        if is_clabsi_suspect:
            verdict = AdjudicationVerdict.CLABSI_SUSPECTED
            is_contamination = False
            action = "Consider catheter removal, paired repeat cultures, and targeted therapy"
        elif score >= 2:
            verdict = AdjudicationVerdict.TRUE_PATHOGEN
            is_contamination = False
            action = "Treat as true bacteremia; optimize targeted antimicrobial therapy"
        elif score <= -2:
            verdict = AdjudicationVerdict.PROBABLE_CONTAMINATION
            is_contamination = True
            action = "Probable contamination; avoid unnecessary vancomycin/antimicrobial escalation"
        else:
            verdict = AdjudicationVerdict.INDETERMINATE
            is_contamination = False
            action = "Indeterminate; evaluate clinical signs, inflammatory markers (CRP/PCT), and repeat culture"

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


# ==============================================================================
# 2. STATISTICAL SURVEILLANCE & WILSON SCORE CONFIDENCE INTERVALS
# ==============================================================================

class ContaminationSurveillanceEngine:
    """Computes contamination rates, Wilson confidence intervals, and unit-level funnel limits."""

    CLSI_TARGET_PCT: float = 3.0  # CLSI standard: < 3.0% contamination rate target

    @classmethod
    def calculate_wilson_ci(cls, contaminated_count: int, total_sets: int, z_score: float = 1.96) -> Tuple[float, float]:
        """
        Calculate 95% Wilson Score Confidence Interval for a proportion:
        p_hat = k/n
        denominator = 1 + z^2 / n
        center = (p_hat + z^2 / (2n)) / denominator
        half_width = z * sqrt((p_hat*(1 - p_hat) + z^2 / (4n)) / n) / denominator
        """
        if total_sets <= 0:
            return 0.0, 0.0

        p = contaminated_count / total_sets
        z2 = z_score * z_score
        n = total_sets

        denom = 1.0 + z2 / n
        center = (p + z2 / (2.0 * n)) / denom
        margin = z_score * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n) / denom

        lower = max(0.0, (center - margin) * 100.0)
        upper = min(100.0, (center + margin) * 100.0)
        return round(lower, 2), round(upper, 2)

    @classmethod
    def analyze_surveillance_data(
        cls,
        culture_records: List[Dict[str, Any]],
        target_pct: float = CLSI_TARGET_PCT,
    ) -> Dict[str, Any]:
        """
        Analyze multi-unit hospital dataset for contamination metrics and SPC funnel control limits.
        """
        total_sets = len(culture_records)
        contaminated_sets = sum(1 for r in culture_records if r.get("is_contamination", False))

        rate_pct = (contaminated_sets / total_sets * 100.0) if total_sets > 0 else 0.0
        ci_lower, ci_upper = cls.calculate_wilson_ci(contaminated_sets, total_sets)

        # Unit breakdown and funnel limits
        unit_buckets: Dict[str, Dict[str, int]] = {}
        for r in culture_records:
            unit = r.get("collection_unit", "Unknown")
            agg = unit_buckets.setdefault(unit, {"total": 0, "contaminated": 0})
            agg["total"] += 1
            if r.get("is_contamination", False):
                agg["contaminated"] += 1

        overall_p = contaminated_sets / max(total_sets, 1)
        funnel_analysis = {}

        for unit_name, counts in unit_buckets.items():
            n_u = counts["total"]
            k_u = counts["contaminated"]
            u_rate = (k_u / n_u * 100.0) if n_u > 0 else 0.0

            # 3-sigma (99.7%) and 2-sigma (95%) upper funnel control limits
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

        meets_standard = rate_pct <= target_pct

        return {
            "total_sets_evaluated": total_sets,
            "contaminated_sets_count": contaminated_sets,
            "overall_contamination_rate_pct": round(rate_pct, 2),
            "wilson_95ci_pct": [ci_lower, ci_upper],
            "clsi_target_threshold_pct": target_pct,
            "meets_clsi_standard": meets_standard,
            "status_summary": "Compliant (Within CLSI Benchmark)" if meets_standard else "NON-COMPLIANT (Exceeds CLSI Benchmark)",
            "unit_funnel_analytics": funnel_analysis,
        }


# ==============================================================================
# 3. ECONOMIC IMPACT MODELING
# ==============================================================================

class EconomicImpactEngine:
    """Calculates excess hospital costs and avoidable antimicrobial expenditures."""

    @classmethod
    def calculate_cost(
        cls,
        contaminated_count: int,
        custom_cost_weights: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """
        Quantify economic burden of contaminated blood culture sets.
        """
        costs = dict(DEFAULT_COST_PER_CONTAMINATION_USD)
        if custom_cost_weights:
            costs.update(custom_cost_weights)

        per_case_total = sum(costs.values())
        total_cost = per_case_total * contaminated_count
        avoidable_if_reduced_by_half = total_cost * 0.50

        return {
            "contaminated_sets_count": contaminated_count,
            "cost_breakdown_per_contamination_usd": costs,
            "total_cost_per_case_usd": round(per_case_total, 2),
            "estimated_annual_excess_cost_usd": round(total_cost, 2),
            "potential_savings_50pct_reduction_usd": round(avoidable_if_reduced_by_half, 2),
        }
