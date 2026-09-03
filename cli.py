"""
Clinical Surveillance & Adjudication CLI for Blood Culture Contamination Tracker
Adheres to CLSI M47-A2, Richter et al. probability algorithm, and CDC NHSN benchmarks.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from blood_culture_tracker import (
    COMMON_CONTAMINANTS,
    DEFAULT_COST_PER_CONTAMINATION_USD,
    TRUE_PATHOGENS,
    BloodCultureSet,
    ContaminationSurveillanceEngine,
    CultureAdjudicationEngine,
    CultureBottle,
    EconomicImpactEngine,
)


def adjudicate_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """
    Parses a single row representing a blood culture draw/set and executes
    multi-criterion adjudication:
    1. Organism pathogenicity (CLSI/Richter criteria: true pathogen vs skin commensal).
    2. Number of positive bottles/sets out of total drawn (concordance).
    3. Time-to-positivity (TTP) and differential TTP (DTTP between line and peripheral).
    4. Phlebotomy vs IV line draw site risk stratification.
    """
    set_id = row.get("set_id") or row.get("culture_id") or row.get("id") or "SET-001"
    patient_id = row.get("patient_id") or row.get("mrn") or "PT-UNKNOWN"
    unit = row.get("collection_unit") or row.get("unit") or "Inpatient"
    organism_raw = row.get("organism") or row.get("name") or row.get("query") or "unknown"
    organism_clean = organism_raw.strip().lower().replace(" ", "_")

    # Positive and total bottles drawn
    try:
        bottles_drawn = int(row.get("bottles_drawn") or row.get("total_bottles") or 2)
    except (ValueError, TypeError):
        bottles_drawn = 2

    try:
        bottles_positive = int(row.get("bottles_positive") or row.get("positive_bottles") or 1)
    except (ValueError, TypeError):
        bottles_positive = 1

    # TTP values
    def parse_float(val: Any) -> Optional[float]:
        if val is None or str(val).strip() == "" or str(val).lower() == "none":
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    ttp_hours = parse_float(row.get("ttp_hours") or row.get("ttp"))
    peripheral_ttp = parse_float(row.get("peripheral_ttp_hours") or row.get("peripheral_ttp"))
    central_ttp = parse_float(row.get("central_ttp_hours") or row.get("central_ttp"))
    draw_site = (row.get("draw_site") or row.get("site_type") or "peripheral").strip().lower()

    bottles = []
    # If explicit central/peripheral TTPs are supplied, instantiate matching bottles
    if peripheral_ttp is not None:
        bottles.append(CultureBottle(
            bottle_id=f"{set_id}-PERIPH",
            site_type="peripheral",
            ttp_hours=peripheral_ttp,
            is_positive=True,
            organism=organism_clean,
        ))
    if central_ttp is not None:
        bottles.append(CultureBottle(
            bottle_id=f"{set_id}-CENTRAL",
            site_type="central_line",
            ttp_hours=central_ttp,
            is_positive=True,
            organism=organism_clean,
        ))

    # Fill remaining bottles according to counts
    needed_pos = max(0, bottles_positive - len(bottles))
    for i in range(needed_pos):
        bottles.append(CultureBottle(
            bottle_id=f"{set_id}-POS-{i+1}",
            site_type="central_line" if "line" in draw_site or "catheter" in draw_site or "picc" in draw_site else "peripheral",
            ttp_hours=ttp_hours,
            is_positive=True,
            organism=organism_clean,
        ))

    needed_neg = max(0, bottles_drawn - len(bottles))
    for i in range(needed_neg):
        bottles.append(CultureBottle(
            bottle_id=f"{set_id}-NEG-{i+1}",
            site_type="peripheral",
            ttp_hours=None,
            is_positive=False,
            organism=None,
        ))

    culture_set = BloodCultureSet(
        set_id=set_id,
        patient_id=patient_id,
        collection_unit=unit,
        bottles=bottles,
    )

    result = CultureAdjudicationEngine.adjudicate_set(culture_set)

    # Calculate Richter-style probability estimate
    # True pathogen: high prior (0.95+); Contaminant species: 0.15-0.80 depending on bottle concordance
    is_pathogen = organism_clean in TRUE_PATHOGENS or any(p in organism_clean for p in TRUE_PATHOGENS)
    is_contaminant = organism_clean in COMMON_CONTAMINANTS or any(c in organism_clean for c in COMMON_CONTAMINANTS)
    concordance = bottles_positive / max(bottles_drawn, 1)

    if is_pathogen:
        prob_true = min(0.99, 0.85 + (0.14 * concordance))
    elif is_contaminant:
        if concordance <= 0.25:
            prob_true = 0.05
        elif concordance <= 0.5:
            prob_true = 0.20 if (ttp_hours and ttp_hours > 30) else 0.35
        elif concordance <= 0.75:
            prob_true = 0.55
        else:
            prob_true = 0.82
    else:
        prob_true = 0.50

    # Draw site risk score: IV/central catheter line draws have higher contamination risk
    # when skin flora is present unless DTTP >= 2.0 hrs confirms CLABSI
    line_draw_flag = "line" in draw_site or "catheter" in draw_site or "picc" in draw_site
    site_risk = "Elevated (Catheter/IV Draw)" if line_draw_flag else "Standard (Venipuncture)"

    return {
        "set_id": set_id,
        "patient_id": patient_id,
        "collection_unit": unit,
        "organism": organism_raw,
        "draw_site": draw_site,
        "site_risk": site_risk,
        "bottles_drawn": bottles_drawn,
        "bottles_positive": bottles_positive,
        "concordance_ratio": round(concordance, 2),
        "ttp_hours": ttp_hours if ttp_hours is not None else "",
        "dttp_hours": result.get("dttp_hours") if result.get("dttp_hours") is not None else "",
        "adjudication_verdict": result["verdict"],
        "is_contamination": result["is_contamination"],
        "confidence_score": result["confidence_score"],
        "true_bacteremia_probability": round(prob_true, 3),
        "clinical_action": result["clinical_recommendation"],
    }


def process_csv(input_path: str, output_path: str) -> List[Dict[str, Any]]:
    """Batch-processes blood culture records from CSV and saves detailed adjudication."""
    with open(input_path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    adjudicated_records: List[Dict[str, Any]] = []
    for r in rows:
        adj = adjudicate_row(r)
        # Merge input fields with adjudication output
        merged = {**r, **adj}
        adjudicated_records.append(merged)

    # Determine CSV output fields
    ordered_fields = list(fieldnames)
    for k in [
        "adjudication_verdict",
        "is_contamination",
        "confidence_score",
        "true_bacteremia_probability",
        "site_risk",
        "concordance_ratio",
        "dttp_hours",
        "clinical_action",
    ]:
        if k not in ordered_fields:
            ordered_fields.append(k)

    with open(output_path, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ordered_fields)
        writer.writeheader()
        writer.writerows(adjudicated_records)

    return adjudicated_records


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="blood-culture-contamination-tracker",
        description="Clinical Microbiology Blood Culture Contamination Tracker (CLSI M47-A2 / CDC NHSN / CAP)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Batch command
    p_batch = subparsers.add_parser("batch", help="Batch process blood culture surveillance CSV records")
    p_batch.add_argument("-i", "--input", required=True, help="Path to input surveillance CSV")
    p_batch.add_argument("-o", "--output", default="results.csv", help="Path to output adjudication CSV")
    p_batch.add_argument("--benchmark-target", type=float, default=3.0, help="Target contamination rate benchmark (default: 3.0%%)")

    # Adjudicate single command
    p_adj = subparsers.add_parser("adjudicate", help="Adjudicate a single blood culture set")
    p_adj.add_argument("--set-id", default="BC-TEST-01", help="Culture set identifier")
    p_adj.add_argument("--organism", required=True, help="Identified bacterial/fungal organism")
    p_adj.add_argument("--bottles-drawn", type=int, default=2, help="Total bottles drawn (default: 2)")
    p_adj.add_argument("--bottles-positive", type=int, default=1, help="Positive bottles count")
    p_adj.add_argument("--ttp", type=float, default=None, help="Time to positivity in hours")
    p_adj.add_argument("--peripheral-ttp", type=float, default=None, help="Peripheral blood TTP (hours)")
    p_adj.add_argument("--central-ttp", type=float, default=None, help="Central line TTP (hours)")
    p_adj.add_argument("--site", default="peripheral", choices=["peripheral", "central_line", "picc"], help="Collection site")

    # Rate surveillance analytics
    p_surv = subparsers.add_parser("surveillance", help="Calculate institutional contamination rate and Wilson 95%% CI")
    p_surv.add_argument("-i", "--input", required=True, help="Path to processed CSV containing 'is_contamination' column")
    p_surv.add_argument("--target", type=float, default=3.0, help="CLSI benchmark threshold %% (default 3.0)")

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "batch":
        records = process_csv(args.input, args.output)
        surv_data = [
            {"is_contamination": r.get("is_contamination", False), "collection_unit": r.get("collection_unit", "Unknown")}
            for r in records
        ]
        stats = ContaminationSurveillanceEngine.analyze_surveillance_data(surv_data, target_pct=args.benchmark_target)
        econ = EconomicImpactEngine.calculate_cost(stats["contaminated_sets_count"])

        print("=" * 80)
        print("  CLSI M47-A2 BLOOD CULTURE CONTAMINATION SURVEILLANCE REPORT")
        print("=" * 80)
        print(f"  Total Sets Evaluated:       {stats['total_sets_evaluated']}")
        print(f"  Contaminated Sets:          {stats['contaminated_sets_count']}")
        print(f"  Institutional Rate:         {stats['overall_contamination_rate_pct']:.2f}% (Target: <= {stats['clsi_target_threshold_pct']:.1f}%)")
        print(f"  Wilson 95% Confidence Int:  [{stats['wilson_95ci_pct'][0]:.2f}%, {stats['wilson_95ci_pct'][1]:.2f}%]")
        print(f"  Benchmark Status:           {stats['status_summary']}")
        print(f"  Estimated Excess Costs:     ${econ['estimated_annual_excess_cost_usd']:,.2f}")
        print(f"  Potential Savings (50% red):${econ['potential_savings_50pct_reduction_usd']:,.2f}")
        print("-" * 80)
        print(f"  Detailed records saved to:   {args.output}")
        print("=" * 80)
        return 0

    if args.command == "adjudicate":
        row_input = {
            "set_id": args.set_id,
            "organism": args.organism,
            "bottles_drawn": args.bottles_drawn,
            "bottles_positive": args.bottles_positive,
            "ttp_hours": args.ttp,
            "peripheral_ttp_hours": args.peripheral_ttp,
            "central_ttp_hours": args.central_ttp,
            "draw_site": args.site,
        }
        res = adjudicate_row(row_input)
        print("=" * 80)
        print("  BLOOD CULTURE ADJUDICATION RESULT")
        print("=" * 80)
        print(f"  Set ID:            {res['set_id']}")
        print(f"  Organism:          {res['organism']}")
        print(f"  Positivity Ratio:  {res['bottles_positive']}/{res['bottles_drawn']} (Concordance: {res['concordance_ratio']})")
        if res['ttp_hours']:
            print(f"  Shortest TTP:      {res['ttp_hours']} hours")
        if res['dttp_hours']:
            print(f"  Differential TTP:  {res['dttp_hours']} hours")
        print(f"  Draw Site Risk:    {res['site_risk']}")
        print(f"  Adjudication:      {res['adjudication_verdict']}")
        print(f"  True Bacteremia P: {res['true_bacteremia_probability']:.1%}")
        print(f"  Contaminant Flag:  {res['is_contamination']}")
        print(f"  Recommendation:    {res['clinical_action']}")
        print("=" * 80)
        return 0

    if args.command == "surveillance":
        with open(args.input, mode="r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        surv_data = []
        for r in rows:
            is_c = str(r.get("is_contamination", "")).strip().lower() in ("true", "1", "yes")
            surv_data.append({"is_contamination": is_c, "collection_unit": r.get("collection_unit", "Unknown")})
        stats = ContaminationSurveillanceEngine.analyze_surveillance_data(surv_data, target_pct=args.target)
        print(json.dumps(stats, indent=2))
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
