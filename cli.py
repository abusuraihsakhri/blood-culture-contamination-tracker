"""Command-line interface for the blood culture contamination tracker."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import Any, Dict, List, Optional

from blood_culture_tracker import (
    COMMON_CONTAMINANTS,
    TRUE_PATHOGENS,
    BloodCultureSet,
    ContaminationSurveillanceEngine,
    CultureAdjudicationEngine,
    CultureBottle,
    EconomicImpactEngine,
    normalize_organism,
)


def _first_present(row: Dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        if key in row and row[key] is not None and str(row[key]).strip() != "":
            return row[key]
    return default


def _parse_float(value: Any) -> Optional[float]:
    if value is None or str(value).strip().lower() in {"", "none", "null"}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric value: {value!r}") from exc


def _parse_int(value: Any, *, default: int, name: str) -> int:
    if value is None or str(value).strip() == "":
        return default
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer.") from exc


def adjudicate_row(row: Dict[str, Any]) -> Dict[str, Any]:
    set_id = str(_first_present(row, "set_id", "culture_id", "id", default="SET-001"))
    patient_id = str(_first_present(row, "patient_id", "mrn", default="PT-UNKNOWN"))
    unit = str(_first_present(row, "collection_unit", "unit", default="Inpatient"))
    organism_raw = str(_first_present(row, "organism", "name", "query", default="unknown"))
    organism_clean = normalize_organism(organism_raw) or "unknown"

    bottles_drawn = _parse_int(
        _first_present(row, "bottles_drawn", "total_bottles"),
        default=2,
        name="bottles_drawn",
    )
    bottles_positive = _parse_int(
        _first_present(row, "bottles_positive", "positive_bottles"),
        default=1,
        name="bottles_positive",
    )
    if bottles_drawn <= 0:
        raise ValueError("bottles_drawn must be greater than zero.")
    if not 0 <= bottles_positive <= bottles_drawn:
        raise ValueError("bottles_positive must satisfy 0 <= bottles_positive <= bottles_drawn.")

    ttp_hours = _parse_float(_first_present(row, "ttp_hours", "ttp"))
    peripheral_ttp = _parse_float(_first_present(row, "peripheral_ttp_hours", "peripheral_ttp"))
    central_ttp = _parse_float(_first_present(row, "central_ttp_hours", "central_ttp"))
    for label, value in (
        ("ttp_hours", ttp_hours),
        ("peripheral_ttp_hours", peripheral_ttp),
        ("central_ttp_hours", central_ttp),
    ):
        if value is not None and value < 0:
            raise ValueError(f"{label} must be non-negative.")

    draw_site = str(_first_present(row, "draw_site", "site_type", default="peripheral")).strip().lower()
    if draw_site not in {"peripheral", "central_line", "picc", "catheter"}:
        raise ValueError("draw_site must be peripheral, central_line, picc, or catheter.")

    explicit_positive_count = int(peripheral_ttp is not None) + int(central_ttp is not None)
    if explicit_positive_count > bottles_positive:
        raise ValueError(
            "Paired TTP inputs imply more positive bottles than bottles_positive reports."
        )

    bottles: List[CultureBottle] = []
    if peripheral_ttp is not None:
        bottles.append(
            CultureBottle(
                bottle_id=f"{set_id}-PERIPH",
                site_type="peripheral",
                ttp_hours=peripheral_ttp,
                is_positive=True,
                organism=organism_clean,
            )
        )
    if central_ttp is not None:
        bottles.append(
            CultureBottle(
                bottle_id=f"{set_id}-CENTRAL",
                site_type="central_line",
                ttp_hours=central_ttp,
                is_positive=True,
                organism=organism_clean,
            )
        )

    positive_site = "central_line" if draw_site in {"central_line", "picc", "catheter"} else "peripheral"
    for i in range(max(0, bottles_positive - len(bottles))):
        bottles.append(
            CultureBottle(
                bottle_id=f"{set_id}-POS-{i + 1}",
                site_type=positive_site,
                ttp_hours=ttp_hours,
                is_positive=True,
                organism=organism_clean,
            )
        )

    for i in range(max(0, bottles_drawn - len(bottles))):
        bottles.append(
            CultureBottle(
                bottle_id=f"{set_id}-NEG-{i + 1}",
                site_type="peripheral",
                is_positive=False,
            )
        )

    culture_set = BloodCultureSet(
        set_id=set_id,
        patient_id=patient_id,
        collection_unit=unit,
        bottles=bottles,
    )
    result = CultureAdjudicationEngine.adjudicate_set(culture_set)

    is_pathogen = organism_clean in TRUE_PATHOGENS
    is_contaminant = organism_clean in COMMON_CONTAMINANTS
    concordance = bottles_positive / bottles_drawn

    # Backward-compatible heuristic probability. It is deliberately labelled as
    # a heuristic and must not be interpreted as a validated clinical model.
    if bottles_positive == 0:
        prob_true = 0.0
    elif is_pathogen:
        prob_true = min(0.99, 0.85 + 0.14 * concordance)
    elif is_contaminant:
        if concordance <= 0.25:
            prob_true = 0.05
        elif concordance <= 0.5:
            prob_true = 0.20 if (ttp_hours is not None and ttp_hours > 30) else 0.35
        elif concordance <= 0.75:
            prob_true = 0.55
        else:
            prob_true = 0.82
    else:
        prob_true = 0.50

    line_draw_flag = draw_site in {"central_line", "picc", "catheter"}
    return {
        "set_id": set_id,
        "patient_id": patient_id,
        "collection_unit": unit,
        "organism": organism_raw,
        "draw_site": draw_site,
        "site_risk": "Higher contamination risk (catheter draw)" if line_draw_flag else "Peripheral venipuncture",
        "bottles_drawn": bottles_drawn,
        "bottles_positive": bottles_positive,
        "concordance_ratio": round(concordance, 3),
        "ttp_hours": "" if ttp_hours is None else ttp_hours,
        "dttp_hours": "" if result.get("dttp_hours") is None else result["dttp_hours"],
        "adjudication_verdict": result["verdict"],
        "is_contamination": result["is_contamination"],
        "confidence_score": result["confidence_score"],
        "true_bacteremia_probability": round(prob_true, 3),
        "probability_note": "Heuristic, not a validated probability model.",
        "clinical_action": result["clinical_recommendation"],
    }


def process_csv(input_path: str, output_path: str) -> List[Dict[str, Any]]:
    with open(input_path, mode="r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        input_fields = list(reader.fieldnames or [])

    adjudicated_records = [{**row, **adjudicate_row(row)} for row in rows]

    ordered_fields: List[str] = list(input_fields)
    for record in adjudicated_records:
        for key in record:
            if key not in ordered_fields:
                ordered_fields.append(key)

    with open(output_path, mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        writer.writerows(adjudicated_records)

    return adjudicated_records


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="blood-culture-contamination-tracker",
        description="Blood-culture contamination and DTTP decision-support utilities",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_batch = subparsers.add_parser("batch", help="Batch process surveillance CSV records")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV")
    p_batch.add_argument("-o", "--output", default="results.csv", help="Output CSV")
    p_batch.add_argument(
        "--benchmark-target",
        type=float,
        default=3.0,
        help="Configured contamination-rate target percentage (default: 3.0)",
    )

    p_adj = subparsers.add_parser("adjudicate", help="Adjudicate one represented culture set")
    p_adj.add_argument("--set-id", default="BC-TEST-01")
    p_adj.add_argument("--organism", required=True)
    p_adj.add_argument("--bottles-drawn", type=int, default=2)
    p_adj.add_argument("--bottles-positive", type=int, default=1)
    p_adj.add_argument("--ttp", type=float, default=None)
    p_adj.add_argument("--peripheral-ttp", type=float, default=None)
    p_adj.add_argument("--central-ttp", type=float, default=None)
    p_adj.add_argument(
        "--site",
        default="peripheral",
        choices=["peripheral", "central_line", "picc", "catheter"],
    )

    p_surv = subparsers.add_parser("surveillance", help="Calculate rate and Wilson 95% CI")
    p_surv.add_argument("-i", "--input", required=True)
    p_surv.add_argument("--target", type=float, default=3.0)

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "batch":
        records = process_csv(args.input, args.output)
        surveillance_rows = [
            {
                "is_contamination": record.get("is_contamination", False),
                "collection_unit": record.get("collection_unit", "Unknown"),
            }
            for record in records
        ]
        stats = ContaminationSurveillanceEngine.analyze_surveillance_data(
            surveillance_rows, target_pct=args.benchmark_target
        )
        econ = EconomicImpactEngine.calculate_cost(stats["contaminated_sets_count"])

        print("=" * 78)
        print("  BLOOD CULTURE CONTAMINATION SURVEILLANCE")
        print("=" * 78)
        print(f"  Total sets:             {stats['total_sets_evaluated']}")
        print(f"  Contaminated sets:      {stats['contaminated_sets_count']}")
        print(
            f"  Contamination rate:     {stats['overall_contamination_rate_pct']:.2f}% "
            f"(configured target <= {stats['configured_target_threshold_pct']:.1f}%)"
        )
        print(
            f"  Wilson 95% CI:          [{stats['wilson_95ci_pct'][0]:.2f}%, "
            f"{stats['wilson_95ci_pct'][1]:.2f}%]"
        )
        print(f"  Target status:          {stats['status_summary']}")
        print(
            f"  Illustrative cost:      ${econ['estimated_annual_excess_cost_usd']:,.2f} "
            "(replace defaults with local values)"
        )
        print(f"  Output:                 {args.output}")
        print("=" * 78)
        return 0

    if args.command == "adjudicate":
        result = adjudicate_row(
            {
                "set_id": args.set_id,
                "organism": args.organism,
                "bottles_drawn": args.bottles_drawn,
                "bottles_positive": args.bottles_positive,
                "ttp_hours": args.ttp,
                "peripheral_ttp_hours": args.peripheral_ttp,
                "central_ttp_hours": args.central_ttp,
                "draw_site": args.site,
            }
        )
        print(json.dumps(result, indent=2))
        return 0

    if args.command == "surveillance":
        with open(args.input, mode="r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        surveillance_rows = [
            {
                "is_contamination": str(row.get("is_contamination", "")).strip().lower()
                in {"true", "1", "yes"},
                "collection_unit": row.get("collection_unit") or "Unknown",
            }
            for row in rows
        ]
        print(
            json.dumps(
                ContaminationSurveillanceEngine.analyze_surveillance_data(
                    surveillance_rows, target_pct=args.target
                ),
                indent=2,
            )
        )
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
