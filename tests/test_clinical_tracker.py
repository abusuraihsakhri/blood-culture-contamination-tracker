from pathlib import Path

import pytest

from blood_culture_tracker import (
    AdjudicationVerdict,
    BloodCultureSet,
    ContaminationSurveillanceEngine,
    CultureAdjudicationEngine,
    CultureBottle,
    EconomicImpactEngine,
    normalize_organism,
)
from cli import adjudicate_row, main, process_csv
from simulator import run_simulation
from ttp_differential import CultureSet, adjudicate as adjudicate_ttp


def _set(*bottles):
    return BloodCultureSet("S", "P", "ICU", list(bottles))


def test_true_pathogen_signal():
    result = CultureAdjudicationEngine.adjudicate_set(
        _set(
            CultureBottle("B1", "peripheral", 9.5, True, "Staphylococcus aureus"),
            CultureBottle("B2", "peripheral", 10.0, True, "staphylococcus_aureus"),
        )
    )
    assert result["verdict"] == AdjudicationVerdict.TRUE_PATHOGEN.value
    assert result["is_contamination"] is False


def test_common_commensal_signal():
    result = CultureAdjudicationEngine.adjudicate_set(
        _set(
            CultureBottle("B1", "peripheral", 44.0, True, "coagulase negative staphylococci"),
            CultureBottle("B2", "peripheral", None, False, None),
            CultureBottle("B3", "peripheral", None, False, None),
            CultureBottle("B4", "peripheral", None, False, None),
        )
    )
    assert result["verdict"] == AdjudicationVerdict.PROBABLE_CONTAMINATION.value
    assert result["is_contamination"] is True


def test_dttp_supports_catheter_source_only_when_central_is_earlier():
    supported = CultureAdjudicationEngine.adjudicate_set(
        _set(
            CultureBottle("BP", "peripheral", 18.0, True, "coagulase_negative_staphylococcus"),
            CultureBottle("BC", "central_line", 14.0, True, "coagulase_negative_staphylococcus"),
        )
    )
    assert supported["dttp_hours"] == 4.0
    assert supported["verdict"] == AdjudicationVerdict.CLABSI_SUSPECTED.value
    assert "CRBSI" in supported["verdict"]

    not_supported = CultureAdjudicationEngine.adjudicate_set(
        _set(
            CultureBottle("BP", "peripheral", 12.0, True, "coagulase_negative_staphylococcus"),
            CultureBottle("BC", "central_line", 16.0, True, "coagulase_negative_staphylococcus"),
        )
    )
    assert not_supported["dttp_hours"] == -4.0
    assert not_supported["verdict"] != AdjudicationVerdict.CLABSI_SUSPECTED.value


def test_dttp_requires_matching_organism():
    result = CultureAdjudicationEngine.adjudicate_set(
        _set(
            CultureBottle("BP", "peripheral", 18.0, True, "staphylococcus_aureus"),
            CultureBottle("BC", "central_line", 10.0, True, "escherichia_coli"),
        )
    )
    assert result["dttp_hours"] is None


def test_zero_positive_culture_is_preserved_and_does_not_crash():
    result = adjudicate_row(
        {
            "set_id": "NEG",
            "organism": "staphylococcus_aureus",
            "bottles_drawn": 2,
            "bottles_positive": 0,
            "draw_site": "peripheral",
        }
    )
    assert result["bottles_positive"] == 0
    assert result["true_bacteremia_probability"] == 0.0
    assert result["is_contamination"] is False
    assert "negative" in result["adjudication_verdict"].lower()


def test_invalid_counts_and_inconsistent_paired_ttp_rejected():
    with pytest.raises(ValueError):
        adjudicate_row({"organism": "e_coli", "bottles_drawn": 2, "bottles_positive": 3})
    with pytest.raises(ValueError):
        adjudicate_row(
            {
                "organism": "e_coli",
                "bottles_drawn": 2,
                "bottles_positive": 1,
                "central_ttp_hours": 12,
                "peripheral_ttp_hours": 14,
            }
        )


def test_normalization_is_conservative():
    assert normalize_organism("Coagulase-negative staphylococci") == "coagulase_negative_staphylococcus"
    assert normalize_organism("not_staphylococcus_aureus_variant") == "not_staphylococcus_aureus_variant"


def test_wilson_interval_and_validation():
    lower, upper = ContaminationSurveillanceEngine.calculate_wilson_ci(3, 100)
    assert 0 < lower < 3 < upper < 10
    with pytest.raises(ValueError):
        ContaminationSurveillanceEngine.calculate_wilson_ci(11, 10)
    with pytest.raises(ValueError):
        ContaminationSurveillanceEngine.analyze_surveillance_data([], target_pct=101)


def test_economic_model_is_explicitly_assumption_based():
    result = EconomicImpactEngine.calculate_cost(5)
    assert result["total_cost_per_case_usd"] == 7100.0
    assert "Illustrative" in result["assumption_note"]
    with pytest.raises(ValueError):
        EconomicImpactEngine.calculate_cost(-1)


def test_compatibility_module_uses_correct_dttp_direction():
    result = adjudicate_ttp(
        CultureSet(
            "X",
            "coagulase_negative_staphylococcus",
            12.0,
            2,
            2,
            peripheral_ttp_hours=10.0,
            central_ttp_hours=14.0,
        )
    )
    assert result["central_peripheral_ttp_differential_h"] == -4.0
    assert "catheter-related" not in result["detail"].lower()


def test_csv_batch_and_cli_smoke(tmp_path):
    sample = Path(__file__).parents[1] / "sample.csv"
    output = tmp_path / "result.csv"
    rows = process_csv(str(sample), str(output))
    assert len(rows) == 15
    assert output.exists()
    assert main(["batch", "-i", str(sample), "-o", str(output)]) == 0
    assert main(["surveillance", "-i", str(output)]) == 0


def test_simulator_exercises_real_path():
    result = run_simulation(50, seed=11, emit=False)
    total = (
        result["contamination"]
        + result["true_bsi_signal"]
        + result["indeterminate"]
        + result["catheter_source"]
    )
    assert total == 50
