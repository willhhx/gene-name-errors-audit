"""SQLite store for record-level results. One row per candidate article, plus
one row per spreadsheet file and one row per flagged cell."""
import sqlite3
from pathlib import Path

import config

SCHEMA = """
CREATE TABLE IF NOT EXISTS articles (
  pmcid INTEGER PRIMARY KEY,
  year INTEGER NOT NULL,
  order_idx INTEGER NOT NULL,          -- position in the seeded random order
  status TEXT NOT NULL,                -- see README: no_excel, scanned, ...
  note TEXT,
  supp_formats TEXT,                   -- e.g. xlsx=2;pdf=1 (images excluded)
  n_excel_files INTEGER DEFAULT 0,
  n_zip_files INTEGER DEFAULT 0,
  article_type TEXT, journal TEXT, publisher TEXT, license TEXT, doi TEXT,
  is_retracted INTEGER,
  n_files_scanned INTEGER DEFAULT 0,
  n_files_unparseable INTEGER DEFAULT 0,
  n_files_too_large INTEGER DEFAULT 0,
  has_gene_list INTEGER DEFAULT 0,
  n_gene_list_files INTEGER DEFAULT 0,
  species TEXT,
  human_only INTEGER DEFAULT 0,
  n_errors INTEGER DEFAULT 0,
  any_error_broad INTEGER DEFAULT 0,   -- dates, serial numbers, scientific notation
  any_error_strict INTEGER DEFAULT 0,  -- dates only
  any_renamed_human INTEGER DEFAULT 0,
  any_non_renamed_human INTEGER DEFAULT 0,
  any_mouse INTEGER DEFAULT 0,
  any_other INTEGER DEFAULT 0,
  scanned_at TEXT
);
CREATE INDEX IF NOT EXISTS ix_articles_year ON articles(year, order_idx);
CREATE TABLE IF NOT EXISTS files (
  pmcid INTEGER, name TEXT, ext TEXT, source TEXT,   -- source: direct | zip
  bytes INTEGER, parse_status TEXT, n_sheets INTEGER,
  n_gene_columns INTEGER, species TEXT, n_cells_checked INTEGER, n_errors INTEGER
);
CREATE INDEX IF NOT EXISTS ix_files_pmcid ON files(pmcid);
CREATE TABLE IF NOT EXISTS errors (
  pmcid INTEGER, file TEXT, sheet TEXT, row INTEGER, col INTEGER,
  species TEXT, raw TEXT, form TEXT, month INTEGER, day INTEGER,
  category TEXT, candidate TEXT
);
CREATE INDEX IF NOT EXISTS ix_errors_pmcid ON errors(pmcid);
CREATE TABLE IF NOT EXISTS run_info (key TEXT PRIMARY KEY, value TEXT);
"""

ARTICLE_COLS = [
    "pmcid", "year", "order_idx", "status", "note", "supp_formats", "n_excel_files",
    "n_zip_files", "article_type", "journal", "publisher", "license", "doi",
    "is_retracted", "n_files_scanned", "n_files_unparseable", "n_files_too_large",
    "has_gene_list", "n_gene_list_files", "species", "human_only", "n_errors",
    "any_error_broad", "any_error_strict", "any_renamed_human",
    "any_non_renamed_human", "any_mouse", "any_other", "scanned_at",
]
FILE_COLS = ["pmcid", "name", "ext", "source", "bytes", "parse_status", "n_sheets",
             "n_gene_columns", "species", "n_cells_checked", "n_errors"]
ERROR_COLS = ["pmcid", "file", "sheet", "row", "col", "species", "raw", "form",
              "month", "day", "category", "candidate"]


def connect(path=None) -> sqlite3.Connection:
    path = Path(path or config.DB_PATH)
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(path), timeout=60)
    con.executescript(SCHEMA)
    return con


def _insert(con, table, cols, rows):
    q = "INSERT INTO %s (%s) VALUES (%s)" % (table, ",".join(cols), ",".join("?" * len(cols)))
    con.executemany(q, [[r.get(c) for c in cols] for r in rows])


def save_result(con, article: dict, files: list, errors: list):
    """Replace everything stored for one article (safe to re-run)."""
    pid = article["pmcid"]
    con.execute("DELETE FROM files WHERE pmcid=?", (pid,))
    con.execute("DELETE FROM errors WHERE pmcid=?", (pid,))
    con.execute("DELETE FROM articles WHERE pmcid=?", (pid,))
    _insert(con, "articles", ARTICLE_COLS, [article])
    _insert(con, "files", FILE_COLS, files)
    _insert(con, "errors", ERROR_COLS, errors)


def set_info(con, key, value):
    con.execute("INSERT OR REPLACE INTO run_info(key, value) VALUES (?,?)", (key, str(value)))
    con.commit()


def eligible_pmcids(con, year, target=None):
    """PMCIDs of the first `target` eligible articles (scanned, with a gene
    list) of a year, in the seeded random order."""
    q = ("SELECT pmcid, order_idx FROM articles WHERE year=? AND status='scanned' "
         "AND has_gene_list=1 ORDER BY order_idx")
    rows = con.execute(q, (year,)).fetchall()
    if target:
        rows = rows[:target]
    return [r[0] for r in rows]


def cutoff_order_idx(con, year, target=None):
    """order_idx of the last article kept for a year, so that the exploratory
    counts over all candidates use the same prefix of the random order."""
    q = ("SELECT order_idx FROM articles WHERE year=? AND status='scanned' "
         "AND has_gene_list=1 ORDER BY order_idx")
    rows = [r[0] for r in con.execute(q, (year,))]
    if not rows:
        return -1
    if target and len(rows) >= target:
        return rows[target - 1]
    return con.execute("SELECT MAX(order_idx) FROM articles WHERE year=?", (year,)).fetchone()[0]
