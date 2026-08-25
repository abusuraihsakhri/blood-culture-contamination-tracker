import blood_culture as m, pathlib
def test_lookup():
    r=m.lookup('creatinine')
    assert 'top_hit' in r and 'score' in r
    pdir=pathlib.Path(__file__).parent
    rows=m.process_csv(str(pdir/'sample.csv'), str(pdir/'tmp_test_out.csv'))
    assert len(rows)>=1
