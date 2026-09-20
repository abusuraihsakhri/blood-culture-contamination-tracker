"""Deterministic stress test for the repository's real adjudication path."""

from __future__ import annotations

import argparse
import random
import time
from typing import Dict

from cli import adjudicate_row


ORGANISMS = [
    "staphylococcus_aureus",
    "escherichia_coli",
    "coagulase_negative_staphylococcus",
    "cutibacterium_acnes",
    "corynebacterium_species",
    "pseudomonas_aeruginosa",
]


def run_simulation(iterations: int = 100, seed: int = 7, emit: bool = True) -> Dict[str, float]:
    if iterations <= 0:
        raise ValueError("iterations must be greater than zero.")

    rng = random.Random(seed)
    counts = {"contamination": 0, "true_bsi_signal": 0, "indeterminate": 0, "catheter_source": 0}
    started = time.perf_counter()

    for index in range(iterations):
        organism = rng.choice(ORGANISMS)
        drawn = rng.choice([2, 4])
        positive = rng.randint(0, drawn)
        site = rng.choice(["peripheral", "central_line"])

        row = {
            "set_id": f"SIM-{index + 1:05d}",
            "organism": organism,
            "bottles_drawn": drawn,
            "bottles_positive": positive,
            "draw_site": site,
        }
        if positive:
            row["ttp_hours"] = round(rng.uniform(6.0, 52.0), 1)

        if positive >= 2 and rng.random() < 0.2:
            central = round(rng.uniform(8.0, 24.0), 1)
            row["central_ttp_hours"] = central
            row["peripheral_ttp_hours"] = round(central + rng.uniform(-3.0, 5.0), 1)

        result = adjudicate_row(row)
        verdict = result["adjudication_verdict"].lower()
        if result["is_contamination"]:
            counts["contamination"] += 1
        elif "catheter-related" in verdict:
            counts["catheter_source"] += 1
        elif "true bloodstream" in verdict:
            counts["true_bsi_signal"] += 1
        else:
            counts["indeterminate"] += 1

    elapsed = time.perf_counter() - started
    summary = {
        "iterations": iterations,
        "seed": seed,
        "elapsed_seconds": round(elapsed, 6),
        "records_per_second": round(iterations / max(elapsed, 1e-9), 1),
        **counts,
    }
    if emit:
        for key, value in summary.items():
            print(f"{key}: {value}")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Stress-test the core adjudication path")
    parser.add_argument("iterations", nargs="?", type=int, default=100)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    run_simulation(args.iterations, seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
