import datetime as dt
import zipfile

import gene_scan
from conftest import HUMAN, MOUSE, make_xls, make_xlsx


def scan(path, ref):
    return gene_scan.scan_workbook(path, ref["sets"], ref["hgnc"])


# ---- single cells ----------------------------------------------------------
def test_classify_cell_forms():
    c = gene_scan.classify_cell
    assert c(dt.datetime(2024, 9, 2)) == ("date", 9, 2)
    assert c("2-Sep") == ("date_string", 9, 2)
    assert c("12-DEC") == ("date_string", 12, 12)
    assert c("5-mar") == ("date_string", 3, 5)
    assert c("43711") == ("serial5", 9, 3)
    assert c(43711) == ("serial5", 9, 3)
    assert c(43711.0) == ("serial5", 9, 3)
    assert c(12345) == ("serial5", None, None)   # too old to be a conversion
    assert c("2.31E+13") == ("scientific", None, None)
    assert c(2.31e19) == ("scientific", None, None)
    assert c("15/03/06")[0] == "date_string"
    # not errors
    assert c("TP53") is None
    assert c("SEPTIN2") is None
    assert c(123) is None
    assert c(1234567) is None       # 7 digits, not a serial
    assert c(0.0043) is None
    assert c(None) is None
    assert c(True) is None
    assert c("") is None
    assert c("MARCHF1") is None


def test_categorise(ref):
    h = ref["hgnc"]
    assert gene_scan.categorise("date", 9, 2, "Hsapiens", h)[0] == "renamed_human"
    assert gene_scan.categorise("date", 3, 1, "Hsapiens", h)[0] == "renamed_human"
    assert gene_scan.categorise("date", 12, 1, "Hsapiens", h)[0] == "renamed_human"
    cat, cand = gene_scan.categorise("date", 10, 4, "Hsapiens", h)
    assert cat == "non_renamed_human" and cand == "OCT4"
    assert gene_scan.categorise("date", 9, 2, "Mmusculus", h)[0] == "mouse"
    assert gene_scan.categorise("date", 9, 2, "Drerio", h)[0] == "other"
    assert gene_scan.categorise("serial5", None, None, "Hsapiens", h)[0] == "other"
    assert gene_scan.categorise("serial5", 3, 1, "Hsapiens", h)[0] == "renamed_human"   # serial of 1-Mar
    assert gene_scan.categorise("date", 7, 20, "Hsapiens", h)[0] == "other"  # no gene


# ---- whole files -----------------------------------------------------------
def test_clean_xlsx_has_gene_list_no_errors(tmp_path, ref, clean_rows):
    p = make_xlsx(tmp_path / "clean.xlsx", ["gene", "value"], clean_rows)
    r = scan(p, ref)
    assert r.parse_status == "ok"
    assert r.has_gene_list and r.errors == []
    assert r.gene_columns[0].species == "Hsapiens"


def test_damaged_xlsx_finds_all_five_types(tmp_path, ref, damaged_rows):
    p = make_xlsx(tmp_path / "bad.xlsx", ["gene", "value"], damaged_rows)
    r = scan(p, ref)
    assert r.has_gene_list
    forms = sorted(e.form for e in r.errors)
    assert forms == ["date", "date", "date", "date_string", "scientific", "serial5"]
    cats = sorted(e.category for e in r.errors)
    assert cats.count("renamed_human") == 3   # SEPT2, MARCH1 and "1-Dec"
    assert cats.count("non_renamed_human") == 1  # OCT4
    assert cats.count("other") == 2
    sept = [e for e in r.errors if e.month == 9][0]
    assert sept.candidate == "SEPT2" and sept.species == "Hsapiens"
    assert sept.col == 1 and sept.row == len(HUMAN) + 2


def test_damaged_xls_same_result_as_xlsx(tmp_path, ref):
    rows = [[g, 1.0] for g in HUMAN] + [
        [dt.datetime(2024, 9, 2), 2.0], [dt.datetime(2024, 3, 1), 2.0], [43711, 2.0],
    ]
    p = make_xls(tmp_path / "bad.xls", ["gene", "value"], rows)
    r = scan(p, ref)
    assert r.parse_status == "ok" and r.has_gene_list
    assert sorted(e.form for e in r.errors) == ["date", "date", "serial5"]
    assert sum(e.category == "renamed_human" for e in r.errors) == 2


def test_mouse_column_species_and_category(tmp_path, ref):
    rows = [[g, 1] for g in MOUSE] + [[dt.datetime(2024, 9, 2), 1]]
    r = scan(make_xlsx(tmp_path / "m.xlsx", ["g", "v"], rows), ref)
    assert r.gene_columns[0].species == "Mmusculus"
    assert r.errors[0].category == "mouse"


def test_dates_outside_gene_column_are_ignored(tmp_path, ref):
    # Column A: gene symbols, column B: genuine dates. Only column A is checked.
    rows = [[g, dt.datetime(2024, 5, i + 1)] for i, g in enumerate(HUMAN)]
    r = scan(make_xlsx(tmp_path / "d.xlsx", ["gene", "date"], rows), ref)
    assert r.has_gene_list and r.errors == []


def test_no_gene_column_means_not_a_gene_list(tmp_path, ref):
    rows = [["sample%d" % i, dt.datetime(2024, 9, 2)] for i in range(30)]
    r = scan(make_xlsx(tmp_path / "n.xlsx", ["s", "d"], rows), ref)
    assert r.parse_status == "ok" and not r.has_gene_list and r.errors == []


def test_fewer_than_five_gene_hits_is_not_gene_column(tmp_path, ref):
    rows = [[g] for g in HUMAN[:4]] + [["foo%d" % i] for i in range(20)]
    r = scan(make_xlsx(tmp_path / "f.xlsx", ["x"], rows), ref)
    assert not r.has_gene_list


def test_first_two_rows_skipped_like_original(tmp_path, ref):
    # header row and first data row are skipped (original deleted 2 lines); row 3 is checked
    rows = [["2-Sep"], ["1-Mar"]] + [[g] for g in HUMAN]
    p = make_xlsx(tmp_path / "h.xlsx", ["43711"], rows)
    r = scan(p, ref)
    assert [e.row for e in r.errors] == [3]


def test_leading_empty_rows_are_ignored(tmp_path, ref):
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws["A4"] = "gene"                      # rows 1-3 empty, header in row 4
    for i, g in enumerate(HUMAN[:6], start=5):
        ws.cell(row=i, column=1, value=g)
    ws.cell(row=11, column=1, value="2-Sep")
    p = tmp_path / "lead.xlsx"
    wb.save(p)
    r = scan(p, ref)
    assert r.has_gene_list and [e.row for e in r.errors] == [11]


def test_gene_column_uses_top_20_rows_only(tmp_path, ref):
    rows = [["x%d" % i] for i in range(25)] + [[g] for g in HUMAN]
    r = scan(make_xlsx(tmp_path / "late.xlsx", ["h"], rows), ref)
    assert not r.has_gene_list          # genes start after row 20, as in the original


def test_second_sheet_is_scanned(tmp_path, ref, clean_rows):
    extra = [["gene"]] + [[g] for g in HUMAN] + [[dt.datetime(2023, 3, 1)]]
    p = make_xlsx(tmp_path / "s.xlsx", ["a", "b"], [["x", 1]], extra_sheet=extra)
    r = scan(p, ref)
    assert r.n_sheets == 2 and len(r.errors) == 1 and r.errors[0].sheet == "second"


def test_garbage_files_are_logged_not_crashing(tmp_path, ref):
    bad = tmp_path / "bad.xlsx"
    bad.write_bytes(b"PK\x03\x04 this is not a real workbook")
    assert scan(bad, ref).parse_status.startswith("unparseable")
    html = tmp_path / "page.xls"
    html.write_text("<html><table><tr><td>TP53</td></tr></table></html>")
    assert scan(html, ref).parse_status == "unparseable:html"
    empty = tmp_path / "e.xls"
    empty.write_bytes(b"")
    assert scan(empty, ref).parse_status == "unparseable:empty"
    old = tmp_path / "trunc.xls"
    old.write_bytes(b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" + b"\x00" * 50)
    assert scan(old, ref).parse_status.startswith("unparseable")


def test_text_file_with_xls_suffix_is_read(tmp_path, ref):
    lines = ["gene\tvalue"] + ["%s\t1" % g for g in HUMAN] + ["2-Sep\t1"]
    p = tmp_path / "t.xls"
    p.write_text("\n".join(lines))
    r = scan(p, ref)
    assert r.parse_status == "ok" and len(r.errors) == 1


def test_too_large_file_is_logged(tmp_path, ref, monkeypatch, clean_rows):
    import config
    p = make_xlsx(tmp_path / "big.xlsx", ["g", "v"], clean_rows)
    monkeypatch.setattr(config, "MAX_FILE_MB", 0)
    assert scan(p, ref).parse_status == "too_large"


def test_zip_members_extracted_flat_and_scanned(tmp_path, ref, zip_with_xlsx):
    out = tmp_path / "out"
    out.mkdir()
    got = list(gene_scan.iter_excel_in_zip(zip_with_xlsx, out))
    assert len(got) == 1 and got[0][0] == "supp/inner.xlsx"
    assert got[0][1].parent == out and got[0][1].name.startswith("zipmember_")
    assert len(scan(got[0][1], ref).errors) == 6


def test_zip_path_traversal_name_is_not_used(tmp_path, ref):
    z = tmp_path / "evil.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("../../evil.xlsx", b"x")
    out = tmp_path / "out"
    out.mkdir()
    got = list(gene_scan.iter_excel_in_zip(z, out))
    assert got[0][1].parent == out
    assert not (tmp_path.parent / "evil.xlsx").exists()


def test_bad_zip_raises_parse_error(tmp_path):
    import pytest
    z = tmp_path / "x.zip"
    z.write_bytes(b"not a zip")
    with pytest.raises(gene_scan.ParseError):
        list(gene_scan.iter_excel_in_zip(z, tmp_path))


def test_original_lists_only_flag(ref):
    import config
    s = gene_scan.load_gene_sets(ref["dir"], ref["hgnc"], original_only=True)
    assert "SEPTIN2" not in s["Hsapiens"]
    s2 = gene_scan.load_gene_sets(ref["dir"], ref["hgnc"])
    assert "SEPTIN2" in s2["Hsapiens"]


def test_macos_resource_fork_stubs_in_zip_are_skipped(tmp_path, ref, damaged_rows):
    x = make_xlsx(tmp_path / "inner.xlsx", ["gene", "value"], damaged_rows)
    z = tmp_path / "mac.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(x, "Source data/fig.xlsx")
        zf.writestr("__MACOSX/Source data/._fig.xlsx", b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X")
        zf.writestr("Source data/._other.xlsx", b"\x00\x05\x16\x07")
    out = tmp_path / "out"
    out.mkdir()
    got = list(gene_scan.iter_excel_in_zip(z, out))
    assert [g[0] for g in got] == ["Source data/fig.xlsx"]


def test_repeated_values_count_once_for_gene_column(tmp_path, ref):
    # A column of the same few gene-like tokens repeated is not a gene list (original grep -c counts distinct)
    rows = [[g] for g in (HUMAN[:3] * 8)] + [["2-Sep"]]
    r = scan(make_xlsx(tmp_path / "rep.xlsx", ["h"], rows), ref)
    assert not r.has_gene_list
    rows = [[g] for g in HUMAN[:5]] * 4 + [["2-Sep"]]
    r = scan(make_xlsx(tmp_path / "five.xlsx", ["h"], rows), ref)
    assert r.has_gene_list and len(r.errors) == 1
