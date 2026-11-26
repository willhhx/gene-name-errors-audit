"""Core scanner: find gene-symbol columns in a spreadsheet and flag cells that
Excel has damaged (dates, five-digit serial numbers, scientific notation).

The detection rules follow the 2021 scanner (Abeysooriya et al., gene_names.sh
in markziemann/GeneNameErrors2020):

* a column is a gene column if more than 4 of its top 20 cells are exact
  matches to the gene symbols of one species;
* in a gene column, a cell is flagged if it looks like a date, a five-digit
  number, or a number in scientific notation.

Differences from the original (all listed in README.md):
* Excel files are read directly in Python (openpyxl, xlrd), so a cell that
  Excel stores as a real date is seen as a date, not as a serial number.
* The first two non-empty rows are skipped, as in the original (header row
  plus the first data row, because the original deleted two lines after export).
* The human gene list also includes current HGNC symbols, so files that use
  the new names (SEPTIN2, MARCHF1) are still recognised as gene lists.
"""
from __future__ import annotations

import csv
import datetime as dt
import io
import math
import re
import time
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from itertools import chain, islice
from pathlib import Path

import config

# ---- the original detection regexes (awk) ----------------------------------
_MONTHS = "JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC"
_MONTHS_TITLE = "|".join(m.title() for m in _MONTHS.split("|"))
_MONTHS_LOWER = _MONTHS.lower()
RE_5DIGIT = re.compile(r"^[0-9]{5}$")
RE_DATE_SLASH = re.compile(r"[0-9][0-9]/[0-9][0-9]/[0-9][0-9]")
RE_DATE_DASH = re.compile(r"[0-9][0-9]-[0-9][0-9]-[0-9][0-9]")
RE_DAY_MON = re.compile(r"[0-9]-(%s|%s|%s)" % (_MONTHS, _MONTHS_TITLE, _MONTHS_LOWER))
RE_SCI = re.compile(r"[0-9]\.[0-9][0-9]E\+[0-9][0-9]")
# for pulling day and month out of a "2-Sep" style string
RE_DAY_MON_STRICT = re.compile(r"^\s*(\d{1,2})-([A-Za-z]{3})\s*$")

MONTH_NUM = {m: i + 1 for i, m in enumerate(_MONTHS.split("|"))}

# How Excel reads gene-like text as a date: month word followed by a number.
_DATE_PREFIXES = {
    "JAN": 1, "JANUARY": 1, "FEB": 2, "FEBRUARY": 2, "MAR": 3, "MARCH": 3,
    "APR": 4, "APRIL": 4, "MAY": 5, "JUN": 6, "JUNE": 6, "JUL": 7, "JULY": 7,
    "AUG": 8, "AUGUST": 8, "SEP": 9, "SEPT": 9, "SEPTEMBER": 9, "OCT": 10,
    "OCTOBER": 10, "NOV": 11, "NOVEMBER": 11, "DEC": 12, "DECEMBER": 12,
}
RE_DATE_SYMBOL = re.compile(
    r"^(%s)(\d{1,2})$" % "|".join(sorted(_DATE_PREFIXES, key=len, reverse=True)),
    re.IGNORECASE,
)

CATEGORIES = ("renamed_human", "non_renamed_human", "mouse", "other")


class ParseError(Exception):
    """File could not be read; the message is the reason that gets logged."""


# ---- reference data --------------------------------------------------------
@dataclass
class HgncInfo:
    approved: set = field(default_factory=set)
    # (month, day) -> list of (symbol, kind, approved_symbol); kind is
    # "previous", "approved" or "alias"
    date_index: dict = field(default_factory=dict)


def load_hgnc(path: Path) -> HgncInfo:
    """Read the HGNC complete set (tab separated)."""
    info = HgncInfo()
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sym = (row.get("symbol") or "").strip()
            if not sym:
                continue
            info.approved.add(sym)
            groups = (
                ([sym], "approved"),
                ((row.get("prev_symbol") or "").split("|"), "previous"),
                ((row.get("alias_symbol") or "").split("|"), "alias"),
            )
            for names, kind in groups:
                for name in names:
                    name = name.strip().strip('"')
                    m = RE_DATE_SYMBOL.match(name)
                    if not m:
                        continue
                    day = int(m.group(2))
                    if not 1 <= day <= 31:
                        continue
                    key = (_DATE_PREFIXES[m.group(1).upper()], day)
                    info.date_index.setdefault(key, []).append((name, kind, sym))
    return info


def load_gene_sets(ref_dir: Path, hgnc: HgncInfo | None = None,
                   original_only: bool = False) -> dict:
    """species -> frozenset of symbols. Files are <Species>_genes in ref_dir."""
    sets = {}
    for path in sorted(Path(ref_dir).glob("*_genes")):
        species = path.name.split("_")[0]
        with open(path, encoding="utf-8", errors="replace") as fh:
            sets[species] = {ln.strip() for ln in fh if ln.strip()}
    if hgnc is not None and not original_only:
        sets.setdefault("Hsapiens", set()).update(hgnc.approved)
    return {k: frozenset(v) for k, v in sets.items()}


# ---- results ---------------------------------------------------------------
@dataclass
class GeneColumn:
    sheet: str
    col: int
    species: str
    hits: int


@dataclass
class CellError:
    sheet: str
    row: int
    col: int
    species: str
    raw: str
    form: str          # date | date_string | serial5 | scientific
    month: int | None
    day: int | None
    category: str      # one of CATEGORIES
    candidate: str     # likely original symbol(s), for date forms


@dataclass
class FileScan:
    parse_status: str = "ok"      # ok | unparseable:<why> | too_large | timeout
    n_sheets: int = 0
    n_cells_checked: int = 0
    gene_columns: list = field(default_factory=list)
    errors: list = field(default_factory=list)

    @property
    def has_gene_list(self) -> bool:
        return bool(self.gene_columns)


# ---- cell classification ---------------------------------------------------
def _serial_to_md(v):
    """Month and day of an Excel date serial, if it falls in a plausible range
    (1990 to 2026, when a gene symbol could have been turned into a date)."""
    if 32874 <= v <= 46387:
        d = dt.date(1899, 12, 30) + dt.timedelta(days=int(v))
        return d.month, d.day
    return None, None


def classify_cell(value):
    """Return (form, month, day) if the cell looks damaged, else None."""
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, dt.datetime) or isinstance(value, dt.date):
        return ("date", value.month, value.day)
    if isinstance(value, dt.time):
        return None
    if isinstance(value, (int, float)):
        v = float(value)
        if math.isnan(v) or math.isinf(v):
            return None
        if v.is_integer() and 10000 <= v <= 99999:
            return ("serial5",) + _serial_to_md(v)
        if abs(v) >= 1e10:
            return ("scientific", None, None)
        return None
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return None
        if RE_5DIGIT.match(s):
            return ("serial5",) + _serial_to_md(int(s))
        if RE_DATE_SLASH.search(s) or RE_DATE_DASH.search(s):
            return ("date_string", None, None)
        if RE_DAY_MON.search(s):
            m = RE_DAY_MON_STRICT.match(s)
            if m and m.group(2).upper() in MONTH_NUM:
                return ("date_string", MONTH_NUM[m.group(2).upper()], int(m.group(1)))
            return ("date_string", None, None)
        if RE_SCI.search(s):
            return ("scientific", None, None)
    return None


def _raw_text(value) -> str:
    if isinstance(value, dt.datetime):
        return value.strftime("%Y-%m-%d")
    if isinstance(value, dt.date):
        return value.isoformat()
    if isinstance(value, float) and value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return str(value)


def categorise(form, month, day, species, hgnc: HgncInfo | None):
    """Return (category, candidate_symbols) for a flagged cell."""
    if form not in ("date", "date_string", "serial5") or month is None:
        return "other", ""
    if species == "Mmusculus":
        return "mouse", ""
    if species != "Hsapiens":
        return "other", ""
    entries = hgnc.date_index.get((month, day), []) if hgnc else []
    if not entries:
        return "other", ""
    prev = [e for e in entries if e[1] == "previous" and e[0] not in hgnc.approved]
    cand = "/".join(dict.fromkeys(e[0] for e in (prev or entries)))[:60]
    if prev:
        return "renamed_human", cand
    return "non_renamed_human", cand


# ---- reading files ---------------------------------------------------------
def sniff_format(path: Path) -> str:
    with open(path, "rb") as fh:
        head = fh.read(2048)
    if head[:2] == b"PK":
        return "xlsx"
    if head[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1":
        return "xls"
    low = head.lower()
    if b"<html" in low or b"<!doctype" in low or b"<table" in low or low.startswith(b"<?xml"):
        return "html"
    if not head:
        return "empty"
    if b"\x00" in head:
        return "binary"
    return "text"


def _rows_openpyxl(ws):
    for row in ws.iter_rows(values_only=True):
        yield row


def _xls_value(cell, datemode, xlrd):
    t = cell.ctype
    if t in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK, xlrd.XL_CELL_ERROR):
        return None
    if t == xlrd.XL_CELL_DATE:
        try:
            return xlrd.xldate.xldate_as_datetime(cell.value, datemode)
        except Exception:  # out-of-range serial: keep as a date-like marker
            try:
                y, mo, d, *_ = xlrd.xldate.xldate_as_tuple(cell.value, datemode)
                return dt.date(max(y, 1), max(mo, 1), max(d, 1))
            except Exception:
                return int(cell.value)
    if t == xlrd.XL_CELL_BOOLEAN:
        return bool(cell.value)
    return cell.value


def _rows_xls(sheet, book, xlrd):
    for r in range(sheet.nrows):
        yield tuple(_xls_value(c, book.datemode, xlrd) for c in sheet.row(r))


def _rows_text(path: Path):
    raw = Path(path).read_bytes()
    text = raw.decode("utf-8", errors="replace")
    first = text.split("\n", 1)[0]
    delim = max(("\t", ",", ";"), key=first.count)
    if first.count(delim) == 0:
        delim = "\t"
    for row in csv.reader(io.StringIO(text), delimiter=delim):
        yield tuple(c if c != "" else None for c in row)


@contextmanager
def open_sheets(path: Path):
    """Yield a list of (sheet_name, row_iterator_factory). Raises ParseError."""
    fmt = sniff_format(path)
    if fmt == "xlsx":
        import openpyxl
        try:
            wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        except Exception as exc:  # bad zip, not a workbook, encrypted...
            raise ParseError("openpyxl:" + type(exc).__name__) from exc
        try:
            yield [(ws.title, (lambda ws=ws: _rows_openpyxl(ws))) for ws in wb.worksheets]
        finally:
            wb.close()
    elif fmt == "xls":
        import xlrd
        try:
            book = xlrd.open_workbook(path, on_demand=True)
        except Exception as exc:
            raise ParseError("xlrd:" + type(exc).__name__) from exc
        try:
            sheets = []
            for i in range(book.nsheets):
                sh = book.sheet_by_index(i)
                sheets.append((sh.name, (lambda sh=sh: _rows_xls(sh, book, xlrd))))
            yield sheets
        finally:
            book.release_resources()
    elif fmt == "text":
        yield [("text", lambda: _rows_text(path))]
    else:
        raise ParseError(fmt)


def _is_empty_row(row):
    return all(v is None or (isinstance(v, str) and not v.strip()) for v in row)


def leading_empty_rows(head_rows):
    """readxl (used by the original scanner) drops empty rows at the top of a
    sheet, so the 'top 20 rows' start at the first non-empty row."""
    n = 0
    for r in head_rows:
        if not _is_empty_row(r):
            break
        n += 1
    return n


def _find_gene_columns(head_rows, gene_sets, sheet_name):
    """Original rule: in the top 20 rows of a column (header row included),
    more than 4 DIFFERENT values exactly match one species' gene symbols."""
    lead = leading_empty_rows(head_rows)
    window = head_rows[lead:lead + config.TOP_N_CELLS]
    ncols = min(max((len(r) for r in window), default=0), 400)
    found = []
    for c in range(ncols):
        # distinct values: the original ran `grep -cxFf top20 genelist`, which counts
        # gene-list entries matched, so a repeated value ("A", "A", "A") counts once
        strings = {r[c].strip() for r in window if c < len(r) and isinstance(r[c], str)}
        if len(strings) < config.MIN_GENE_HITS:
            continue
        best_species, best_hits = None, 0
        for species, gset in gene_sets.items():
            hits = len(strings & gset)
            if hits > best_hits:
                best_species, best_hits = species, hits
        if best_species and best_hits >= config.MIN_GENE_HITS:
            found.append(GeneColumn(sheet_name, c, best_species, best_hits))
    return found


def scan_workbook(path, gene_sets: dict, hgnc: HgncInfo | None,
                  deadline_s: float = 120.0) -> FileScan:
    """Scan one spreadsheet file. Never raises; problems go in parse_status."""
    path = Path(path)
    scan = FileScan()
    t0 = time.monotonic()
    try:
        if path.stat().st_size > config.MAX_FILE_MB * 1024 * 1024:
            scan.parse_status = "too_large"
            return scan
        with open_sheets(path) as sheets:
            for name, factory in sheets:
                scan.n_sheets += 1
                rows = factory()
                head = list(islice(rows, config.HEAD_ROWS))
                cols = _find_gene_columns(head, gene_sets, name)
                if not cols:
                    continue
                scan.gene_columns.extend(cols)
                skip_upto = leading_empty_rows(head) + config.HEADER_ROWS_SKIPPED
                n = 0
                for r_idx, row in enumerate(chain(head, rows), start=1):
                    n += 1
                    if n > config.MAX_ROWS_PER_SHEET:
                        break
                    if n % 2000 == 0 and time.monotonic() - t0 > deadline_s:
                        scan.parse_status = "timeout"
                        return scan
                    if r_idx <= skip_upto:
                        continue
                    for gc in cols:
                        if gc.col >= len(row):
                            continue
                        value = row[gc.col]
                        hit = classify_cell(value)
                        scan.n_cells_checked += 1
                        if hit is None:
                            continue
                        form, month, day = hit
                        cat, cand = categorise(form, month, day, gc.species, hgnc)
                        scan.errors.append(CellError(
                            sheet=name, row=r_idx, col=gc.col + 1, species=gc.species,
                            raw=_raw_text(value)[:40], form=form, month=month, day=day,
                            category=cat, candidate=cand))
    except ParseError as exc:
        scan.parse_status = "unparseable:" + str(exc)
    except Exception as exc:  # corrupt content found while streaming rows
        scan.parse_status = "unparseable:" + type(exc).__name__
    return scan


def iter_excel_in_zip(zip_path: Path, out_dir: Path):
    """Yield (member_name, extracted_path) for Excel files inside a zip.

    Members are streamed to out_dir under a flat generated name (never the name
    stored in the zip), and size limits stop zip bombs.
    """
    total = 0
    n = 0
    try:
        zf = zipfile.ZipFile(zip_path)
    except Exception as exc:
        raise ParseError("zip:" + type(exc).__name__) from exc
    with zf:
        for info in zf.infolist():
            low = info.filename.lower()
            if info.is_dir() or not low.endswith((".xls", ".xlsx")):
                continue
            # macOS zips carry "__MACOSX/._name.xlsx" resource-fork stubs, not spreadsheets
            if low.startswith("__macosx/") or low.rsplit("/", 1)[-1].startswith("._"):
                continue
            if n >= config.MAX_ZIP_MEMBERS:
                break
            total += info.file_size
            if total > config.MAX_ZIP_UNCOMPRESSED_MB * 1024 * 1024:
                break
            if info.file_size > config.MAX_FILE_MB * 1024 * 1024:
                yield info.filename, None   # caller logs too_large
                continue
            dest = Path(out_dir) / ("zipmember_%d%s" % (n, ".xlsx" if low.endswith("x") else ".xls"))
            with zf.open(info) as src, open(dest, "wb") as dst:
                dst.write(src.read())
            n += 1
            yield info.filename, dest
