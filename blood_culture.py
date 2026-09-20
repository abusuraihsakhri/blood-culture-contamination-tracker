#!/usr/bin/env python3
"""Small compatibility CLI for organism lookup and CSV annotation."""

import argparse
import csv
import sys

from blood_culture_tracker import COMMON_CONTAMINANTS, TRUE_PATHOGENS, normalize_organism


def lookup(query, extra=None):
    key = normalize_organism(query)
    if key in TRUE_PATHOGENS:
        classification = "configured pathogen"
        score = 10
    elif key in COMMON_CONTAMINANTS:
        classification = "configured common commensal"
        score = 10
    else:
        classification = "unclassified"
        score = 0

    label = f"{key or 'unknown'} ({classification})"
    return {
        "query": query,
        "top_hit": label,
        "score": score,
        "classification": classification,
        "all": [(score, label)],
    }


def process_csv(inp, out):
    with open(inp, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    lowered = {name.lower(): name for name in fieldnames}
    query_column = next(
        (
            lowered[name]
            for name in ("organism", "query", "name", "test", "code")
            if name in lowered
        ),
        fieldnames[0] if fieldnames else None,
    )
    if query_column is None:
        raise ValueError("Input CSV has no header columns.")

    results = []
    for row in rows:
        result = lookup(row.get(query_column, ""), row)
        results.append(
            {
                **row,
                "top_hit": result["top_hit"],
                "lookup_score": result["score"],
                "classification": result["classification"],
            }
        )

    output_fields = list(fieldnames)
    for name in ("top_hit", "lookup_score", "classification"):
        if name not in output_fields:
            output_fields.append(name)

    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(results)
    return results


def build_parser():
    parser = argparse.ArgumentParser(
        prog="blood_culture",
        description="Blood-culture organism list lookup",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    single = sub.add_parser("single")
    single.add_argument("query", nargs="?", default="staphylococcus_aureus")
    single.add_argument("--query", dest="q2")
    batch = sub.add_parser("batch")
    batch.add_argument("--input", required=True)
    batch.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.cmd == "single":
        print(lookup(args.q2 or args.query))
        return 0
    if args.cmd == "batch":
        result = process_csv(args.input, args.output)
        print(f"Processed {len(result)} rows -> {args.output}")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
