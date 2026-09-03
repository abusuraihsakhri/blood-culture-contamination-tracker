import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from blood_culture_tracker import (
    AdjudicationVerdict,
    BloodCultureSet,
    CultureBottle,
    CultureAdjudicationEngine,
    ContaminationSurveillanceEngine,
    EconomicImpactEngine,
)
from ttp_differential import (
    CultureSet,
    adjudicate as adjudicate_ttp,
    contamination_rate_stats,
    economic_impact,
)
from cli import main, adjudicate_row, process_csv


def test_true_pathogen_adjudication():
    # S. aureus rapid growth (TTP 9.5h), 2/2 bottles positive -> true pathogen
    b1 = CultureBottle(bottle_id="B1", site_type="peripheral", ttp_hours=9.5, is_positive=True, organism="staphylococcus_aureus")
    b2 = CultureBottle(bottle_id="B2", site_type="peripheral", ttp_hours=10.0, is_positive=True, organism="staphylococcus_aureus")
    cset = BloodCultureSet(set_id="S1", patient_id="P1", collection_unit="ED", bottles=[b1, b2])
    res = CultureAdjudicationEngine.adjudicate_set(cset)
    assert res["verdict"] == AdjudicationVerdict.TRUE_PATHOGEN.value
    assert res["is_contamination"] is False
    assert res["confidence_score"] >= 2


def test_skin_contaminant_adjudication():
    # CoNS delayed growth (TTP 44h), 1/4 bottles positive -> probable contamination
    bottles = [
        CultureBottle(bottle_id="B1", site_type="peripheral", ttp_hours=44.0, is_positive=True, organism="coagulase_negative_staphylococcus"),
        CultureBottle(bottle_id="B2", site_type="peripheral", ttp_hours=None, is_positive=False),
        CultureBottle(bottle_id="B3", site_type="peripheral", ttp_hours=None, is_positive=False),
        CultureBottle(bottle_id="B4", site_type="peripheral", ttp_hours=None, is_positive=False),
    ]
    cset = BloodCultureSet(set_id="S2", patient_id="P2", collection_unit="Ward_3B", bottles=bottles)
    res = CultureAdjudicationEngine.adjudicate_set(cset)
    assert res["verdict"] == AdjudicationVerdict.PROBABLE_CONTAMINATION.value
    assert res["is_contamination"] is True
    assert res["confidence_score"] <= -2


def test_clabsi_dttp_criterion():
    # Central line turns positive >= 2.0 hours before peripheral blood
    b_periph = CultureBottle(bottle_id="BP", site_type="peripheral", ttp_hours=18.0, is_positive=True, organism="coagulase_negative_staphylococcus")
    b_central = CultureBottle(bottle_id="BC", site_type="central_line", ttp_hours=14.0, is_positive=True, organism="coagulase_negative_staphylococcus")
    cset = BloodCultureSet(set_id="S3", patient_id="P3", collection_unit="ICU", bottles=[b_periph, b_central])
    res = CultureAdjudicationEngine.adjudicate_set(cset)
    assert res["verdict"] == AdjudicationVerdict.CLABSI_SUSPECTED.value
    assert res["dttp_hours"] == 4.0
    assert res["is_contamination"] is False


def test_wilson_ci_and_surveillance():
    # 3 contaminated sets out of 100 sets = 3%
    lower, upper = ContaminationSurveillanceEngine.calculate_wilson_ci(3, 100)
    assert 0.0 < lower < 3.0
    assert 3.0 < upper < 10.0

    # Surveillance analysis with target <= 3.0%
    surv_records = [{"is_contamination": False, "collection_unit": "ED"}] * 98 + [{"is_contamination": True, "collection_unit": "ED"}] * 2
    res = ContaminationSurveillanceEngine.analyze_surveillance_data(surv_records, target_pct=3.0)
    assert res["meets_clsi_standard"] is True
    assert res["overall_contamination_rate_pct"] == 2.0


def test_economic_impact():
    econ = EconomicImpactEngine.calculate_cost(5)
    assert econ["contaminated_sets_count"] == 5
    assert econ["total_cost_per_case_usd"] == 7100.0
    assert econ["estimated_annual_excess_cost_usd"] == 35500.0
    assert econ["potential_savings_50pct_reduction_usd"] == 17750.0


def test_ttp_differential_module():
    cs = CultureSet(
        set_id="BC01",
        organism="staphylococcus_aureus",
        ttp_hours=8.0,
        bottles_drawn=2,
        bottles_positive=2,
    )
    res = adjudicate_ttp(cs)
    assert res["verdict"] == "likely_true_bacteremia"

    econ = economic_impact(2)
    assert econ["contaminated_sets"] == 2
    assert econ["estimated_total_cost_usd"] > 0


def test_cli_batch_and_adjudicate_commands(tmp_path):
    # Test adjudicate_row logic
    row = {
        "set_id": "TEST-ROW",
        "organism": "escherichia_coli",
        "bottles_drawn": 2,
        "bottles_positive": 2,
        "ttp_hours": 8.0,
    }
    adj = adjudicate_row(row)
    assert adj["adjudication_verdict"] == AdjudicationVerdict.TRUE_PATHOGEN.value
    assert adj["true_bacteremia_probability"] > 0.90

    # Test batch execution via cli.main
    root = Path(__file__).parent.parent
    sample_csv = root / "sample.csv"
    out_csv = tmp_path / "test_out.csv"
    ret = main(["batch", "-i", str(sample_csv), "-o", str(out_csv)])
    assert ret == 0
    assert out_csv.exists()

    # Test single adjudicate CLI command
    ret_adj = main(["adjudicate", "--organism", "staphylococcus_aureus", "--bottles-positive", "2", "--ttp", "9.0"])
    assert ret_adj == 0

    # Test surveillance CLI command
    ret_surv = main(["surveillance", "-i", str(out_csv)])
    assert ret_surv == 0
