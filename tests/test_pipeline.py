"""End-to-end test of sampling + scanning against a fake local 'bucket'."""
import json
import zipfile

import pytest

import config
import pipeline
import store
from conftest import HUMAN, make_xls, make_xlsx

import datetime as dt


def build_bucket(root, ref_dir, monkeypatch):
    """Articles:
    101 research article, xlsx with 2 damaged cells
    102 research article, xls, clean gene list
    103 no supplementary Excel (pdf only)
    104 review article with damaged xlsx (excluded by type)
    105 zip containing a damaged xlsx
    106 xlsx without any gene column
    107 corrupt xlsx
    108 not open access
    109 research article, xlsx, clean gene list
    """
    root.mkdir()
    (root / "metadata").mkdir()

    def add(pmcid, article_type="research-article", files=None, oa=True):
        vdir = root / ("PMC%d.1" % pmcid)
        vdir.mkdir()
        names = []
        for name, writer in (files or {}).items():
            writer(vdir / name)
            names.append(name)
        (vdir / ("PMC%d.1.xml" % pmcid)).write_text(
            '<article xmlns:xlink="x" article-type="%s"><front><journal-meta>'
            "<journal-title>Test Journal</journal-title><publisher><publisher-name>Pub"
            "</publisher-name></publisher></journal-meta></front></article>" % article_type)
        meta = {"pmcid": "PMC%d" % pmcid, "version": 1, "doi": "10/x%d" % pmcid,
                "is_pmc_openaccess": oa, "is_retracted": False, "license_code": "CC BY",
                "media_urls": ["s3://pmc-oa-opendata/PMC%d.1/%s?md5=00" % (pmcid, n) for n in names]}
        (root / "metadata" / ("PMC%d.1.json" % pmcid)).write_text(json.dumps(meta))

    bad = [[g, 1] for g in HUMAN] + [[dt.datetime(2024, 9, 2), 1], [43711, 1]]
    ok = [[g, 1] for g in HUMAN]

    def zipper(path):
        tmp = path.parent / "tmp_inner.xlsx"
        make_xlsx(tmp, ["g", "v"], bad)
        with zipfile.ZipFile(path, "w") as zf:
            zf.write(tmp, "d/inner.xlsx")
        tmp.unlink()

    def corrupt(path):
        path.write_bytes(b"PK\x03\x04 broken")

    add(101, files={"s1.xlsx": lambda p: make_xlsx(p, ["g", "v"], bad), "fig.jpg": lambda p: p.write_bytes(b"x")})
    add(102, files={"s1.xls": lambda p: make_xls(p, ["g", "v"], ok)})
    add(103, files={"s1.pdf": lambda p: p.write_bytes(b"%PDF")})
    add(104, article_type="review-article", files={"s1.xlsx": lambda p: make_xlsx(p, ["g", "v"], bad)})
    add(105, files={"s1.zip": zipper})
    add(106, files={"s1.xlsx": lambda p: make_xlsx(p, ["a", "b"], [["x%d" % i, i] for i in range(30)])})
    add(107, files={"s1.xlsx": corrupt})
    add(108, oa=False, files={"s1.xlsx": lambda p: make_xlsx(p, ["g", "v"], bad)})
    add(109, files={"s1.xlsx": lambda p: make_xlsx(p, ["g", "v"], ok)})
    return [101, 102, 103, 104, 105, 106, 107, 108, 109, 110]  # 110 has no metadata


@pytest.fixture
def setup(tmp_path, ref, monkeypatch):
    monkeypatch.setattr(config, "REF", ref["dir"])
    (ref["dir"] / "hgnc_complete_set.txt").write_text(ref["hgnc_path"].read_text())
    bucket = tmp_path / "bucket"
    ids = build_bucket(bucket, ref["dir"], monkeypatch)
    spec = {"source_dir": str(bucket), "original_only": False, "threads": 4}
    con = store.connect(tmp_path / "t.sqlite")
    return ids, spec, con


def statuses(con):
    return dict(con.execute("SELECT pmcid, status FROM articles"))


def test_full_pass_statuses_and_flags(setup):
    ids, spec, con = setup
    pipeline.run_year(con, 2024, ids, spec, target=999, procs=1, log=lambda *_: None)
    st = statuses(con)
    assert st[101] == "scanned" and st[102] == "scanned" and st[105] == "scanned"
    assert st[103] == "no_excel"
    assert st[104] == "excluded_type"
    assert st[106] == "scanned" and st[107] == "scanned" and st[109] == "scanned"
    assert st[108] == "not_open_access"
    assert st[110] == "not_found"

    def art(pid):
        return con.execute("SELECT has_gene_list, any_error_broad, any_error_strict, "
                           "any_renamed_human, n_errors, human_only, n_files_unparseable "
                           "FROM articles WHERE pmcid=?", (pid,)).fetchone()
    assert art(101) == (1, 1, 1, 1, 2, 1, 0)
    assert art(102) == (1, 0, 0, 0, 0, 1, 0)
    assert art(105) == (1, 1, 1, 1, 2, 1, 0)     # error found inside the zip
    assert art(106)[0] == 0                       # no gene column: not eligible
    assert art(107)[0] == 0 and art(107)[6] == 1  # corrupt file logged, not dropped
    assert art(109) == (1, 0, 0, 0, 0, 1, 0)
    assert pipeline.eligible_orders(con, 2024) != []
    n_err = con.execute("SELECT COUNT(*) FROM errors").fetchone()[0]
    assert n_err == 4


def test_supp_formats_exclude_images(setup):
    ids, spec, con = setup
    pipeline.run_year(con, 2024, ids, spec, target=999, procs=1, log=lambda *_: None)
    assert con.execute("SELECT supp_formats FROM articles WHERE pmcid=101").fetchone()[0] == "xlsx=1"
    assert con.execute("SELECT supp_formats FROM articles WHERE pmcid=103").fetchone()[0] == "pdf=1"


def test_order_is_seeded_and_stable():
    a = pipeline.candidate_order(2024, [5, 3, 9, 1, 7, 100, 42])
    b = pipeline.candidate_order(2024, [42, 100, 7, 1, 9, 3, 5, 5])
    assert a == b and sorted(a) == [1, 3, 5, 7, 9, 42, 100]
    assert a != pipeline.candidate_order(2025, [5, 3, 9, 1, 7, 100, 42])


def test_stops_at_target_and_resume_does_not_redo(setup):
    ids, spec, con = setup
    log = []
    n = pipeline.run_year(con, 2024, ids, spec, target=2, procs=1, chunk=2, log=log.append)
    assert n >= 2
    done_first = con.execute("SELECT COUNT(*) FROM articles").fetchone()[0]
    assert done_first < len(ids)                  # stopped early
    stamp = con.execute("SELECT scanned_at FROM articles ORDER BY order_idx LIMIT 1").fetchone()
    pipeline.run_year(con, 2024, ids, spec, target=2, procs=1, chunk=2, log=lambda *_: None)
    assert con.execute("SELECT COUNT(*) FROM articles").fetchone()[0] == done_first  # nothing new
    pipeline.run_year(con, 2024, ids, spec, target=999, procs=1, log=lambda *_: None)
    assert con.execute("SELECT COUNT(*) FROM articles").fetchone()[0] == len(ids)
    assert con.execute("SELECT scanned_at FROM articles ORDER BY order_idx LIMIT 1").fetchone() == stamp


def test_fetch_errors_are_retried(setup, monkeypatch):
    ids, spec, con = setup
    import pmc_access
    real = pmc_access.LocalSource.metadata
    calls = {"n": 0}

    def flaky(self, pmcid):
        if pmcid == 101 and calls["n"] == 0:
            calls["n"] += 1
            raise RuntimeError("network down")
        return real(self, pmcid)
    monkeypatch.setattr(pmc_access.LocalSource, "metadata", flaky)
    pipeline.run_year(con, 2024, ids, spec, target=999, procs=1, log=lambda *_: None)
    assert statuses(con)[101] == "fetch_error"
    pipeline.run_year(con, 2024, ids, spec, target=999, procs=1, log=lambda *_: None)
    assert statuses(con)[101] == "scanned"
    assert con.execute("SELECT COUNT(*) FROM articles WHERE pmcid=101").fetchone()[0] == 1
