"""Shared test fixtures: tiny gene lists, a tiny HGNC file, and spreadsheets
built on the fly (xlsx with openpyxl, xls with xlwt)."""
import datetime as dt
import sys
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import gene_scan  # noqa: E402

HUMAN = ["TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PTEN", "GAPDH", "ACTB",
         "TNF", "IL6", "CDK2", "AKT1", "STAT3", "NOTCH1", "VEGFA", "SOX2"]
MOUSE = ["Trp53", "Brca1", "Egfr", "Myc", "Kras", "Pten", "Gapdh", "Actb",
         "Tnf", "Il6", "Cdk2", "Akt1"]

HGNC_TSV = (
    "hgnc_id\tsymbol\tname\tprev_symbol\talias_symbol\tstatus\n"
    "HGNC:1\tSEPTIN2\tseptin 2\tSEPT2\t\tApproved\n"
    "HGNC:2\tMARCHF1\tmembrane associated ring-CH-type finger 1\tMARCH1\t\tApproved\n"
    "HGNC:3\tDELEC1\tdeleted in esophageal cancer 1\tDEC1\t\tApproved\n"
    "HGNC:4\tPOU5F1\tPOU class 5 homeobox 1\t\tOCT4\tApproved\n"
    "HGNC:5\tTP53\ttumor protein p53\t\t\tApproved\n"
)


@pytest.fixture(scope="session")
def ref(tmp_path_factory):
    d = tmp_path_factory.mktemp("ref")
    (d / "Hsapiens_genes").write_text("\n".join(HUMAN) + "\n")
    (d / "Mmusculus_genes").write_text("\n".join(MOUSE) + "\n")
    hgnc = d / "hgnc.tsv"
    hgnc.write_text(HGNC_TSV)
    info = gene_scan.load_hgnc(hgnc)
    sets = gene_scan.load_gene_sets(d, info)
    return {"dir": d, "hgnc": info, "sets": sets, "hgnc_path": hgnc}


def make_xlsx(path, header, rows, extra_sheet=None):
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(header)
    for r in rows:
        ws.append(r)
    if extra_sheet:
        ws2 = wb.create_sheet("second")
        for r in extra_sheet:
            ws2.append(r)
    wb.save(path)
    return path


def make_xls(path, header, rows):
    import xlwt
    wb = xlwt.Workbook()
    ws = wb.add_sheet("Sheet1")
    date_style = xlwt.easyxf(num_format_str="D-MMM")
    for c, h in enumerate(header):
        ws.write(0, c, h)
    for r, row in enumerate(rows, start=1):
        for c, v in enumerate(row):
            if isinstance(v, (dt.datetime, dt.date)):
                ws.write(r, c, v, date_style)
            elif v is not None:
                ws.write(r, c, v)
    wb.save(path)
    return path


@pytest.fixture
def clean_rows():
    return [[g, i * 1.5] for i, g in enumerate(HUMAN)]


@pytest.fixture
def damaged_rows():
    """16 real genes plus 5 damaged cells in the same column."""
    rows = [[g, 1.0] for g in HUMAN]
    rows += [
        [dt.datetime(2024, 9, 2), 2.0],      # SEPT2 -> 2-Sep (renamed gene)
        [dt.datetime(2024, 3, 1), 2.0],      # MARCH1 -> 1-Mar (renamed gene)
        [dt.datetime(2024, 10, 4), 2.0],     # OCT4 -> 4-Oct (not renamed)
        ["1-Dec", 2.0],                      # text form of DEC1
        [43711, 2.0],                        # five-digit serial number
        [2.31e19, 2.0],                      # RIKEN-style id turned into 2.31E+19
    ]
    return rows


@pytest.fixture
def zip_with_xlsx(tmp_path, damaged_rows):
    x = make_xlsx(tmp_path / "inner.xlsx", ["gene", "value"], damaged_rows)
    z = tmp_path / "bundle.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(x, "supp/inner.xlsx")
        zf.writestr("readme.txt", "hello")
    return z
