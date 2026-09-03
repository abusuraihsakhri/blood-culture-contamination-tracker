import blood_culture as m, pathlib

def test_lookup():
    r = m.lookup('staphylococcus_aureus')
    assert 'top_hit' in r and 'score' in r
    root = pathlib.Path(__file__).parent.parent
    out_file = root / 'tmp_test_out.csv'
    try:
        rows = m.process_csv(str(root / 'sample.csv'), str(out_file))
        assert len(rows) >= 1
    finally:
        if out_file.exists():
            out_file.unlink()

