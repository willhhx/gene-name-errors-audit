# Excel gene-name errors after the fixes (2021 to 2025)

Code for a cross-sectional audit of supplementary Excel gene lists. Author: Haoxi Hu.
Question: among open-access PMC papers from 2019 to 2025 with Excel gene lists, how many contain
gene names damaged by Excel, and did that share fall after the 2020 gene renaming and the
October 2023 Excel switch?

## Quick start (Mac)

```bash
cd "Python Intro/excel_gene_errors"
pip install -r requirements-dev.txt        # or: conda env create -f environment.yml
python -m pytest tests -q                  # 30 tests, a few seconds
python run_all.py --pilot --procs 4        # about 10 minutes, tiny real-data check
python run_all.py --procs 6                # the full study, hours, safe to stop and restart
```

Optional: `export NCBI_API_KEY=...` makes the candidate search faster (not required).
Needs internet. Uses about a few GB of disk: spreadsheets are deleted right after scanning.

## What each step does

| Step | Script | What it does |
|---|---|---|
| 0 | `00_setup_reference.py` | Downloads the 2021 study's gene lists and raw results, HGNC symbols. Rebuilds the published 2014 to 2020 yearly rates (they reproduce 3,436 of 11,117 = 30.9% exactly). |
| 1 | `01_build_candidates.py` | PMC search `genom*[Title/Abstract] AND open_access[Filter]`, per year (see note below on why not `[Abstract]`). |
| 2 | `02_scan_articles.py` | Seeded random order per year. Downloads Excel files (also inside zips), scans them, stores everything in `data/scan.sqlite`. Stops at 1,500 eligible articles per year. |
| 3 | `03_validation.py` | `make` picks 200 flagged + 200 unflagged files, `fetch` downloads them, files are labelled by hand, `score` gives the positive predictive value and miss rate. |
| 5 | `05_calibrate_original.py` | Rescans a random 250 per year of the 2021 study's own 2019-2020 gene-list articles and measures agreement (kappa). |
| 4 | `04_analyze.py` | Yearly rates (Wilson CI), logistic trend, before/after, error types, sensitivity runs, calibration, Figures 1 and 2, Table 1. Writes `output/`. |

`run_all.py` chains them. Everything is in `config.py` (years, seed, targets, size caps).

## Detection rules (same as the 2021 scanner)

A column is a gene column when more than 4 of its top 20 cells match one species' gene symbols.
In a gene column a cell is flagged if it is a date, a five-digit number, or in scientific notation.
Errors are typed as: renamed human gene (SEPT, MARCH, DEC1 families, from HGNC previous symbols),
other human symbol, mouse symbol, other.

## Search query note

The 2021 study searched `genom*[Abstract]`. In today's PMC search that field behaves like All Fields
(139,408 hits for 2019, ~5 times the original list). `genom*[Title/Abstract]` gives 30,468 and contains
96.6% of the original 2019 list, so it is used.

## Calibration result

Head-to-head on 496 random 2019-2020 articles from the original study: gene list found in 98.8%,
agreement with their confirmed error calls 99.4%, Cohen's kappa 0.99.


## Differences from the 2021 scanner

1. Spreadsheets are read directly (openpyxl, xlrd), so real date cells are seen as dates. The original used R text export.
2. Gene columns and skipped rows follow the original exactly (top 20 rows; more than 4 distinct gene symbols; first 2 rows skipped). Large numbers stored as numbers (e.g. RIKEN IDs) are flagged as scientific notation; the original regex missed some of these.
3. Human gene list = 2015 Ensembl list plus current HGNC symbols, so files with new names (SEPTIN2) are still recognised. Run `python run_all.py --with-original-lists` for the original-lists-only sensitivity run.
4. Five-digit serials are also turned back into a date to type them (1990 to 2026 only).
5. Files over 50 MB are logged as `too_large`, not silently dropped. Corrupt files are logged as `unparseable:<reason>`.
6. Year is the PMC publication date from the search (`[pdat]`).

## Files

`05_calibrate_original.py` calibration · `gene_scan.py` scanner · `pmc_access.py` PMC search and download · `pipeline.py` sampling loop ·
`store.py` SQLite · `tests/`
