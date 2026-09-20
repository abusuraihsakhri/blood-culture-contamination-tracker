from pathlib import Path

import blood_culture as module


def test_lookup_uses_blood_culture_organism_lists():
    pathogen = module.lookup("Staphylococcus aureus")
    contaminant = module.lookup("coagulase-negative staphylococci")
    unknown = module.lookup("example organism")
    assert pathogen["classification"] == "configured pathogen"
    assert contaminant["classification"] == "configured common commensal"
    assert unknown["classification"] == "unclassified"


def test_process_csv_prefers_organism_column(tmp_path):
    root = Path(__file__).parents[1]
    output = tmp_path / "lookup.csv"
    rows = module.process_csv(str(root / "sample.csv"), str(output))
    assert len(rows) == 15
    assert rows[0]["classification"] == "configured common commensal"
    assert rows[1]["classification"] == "configured pathogen"
