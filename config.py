"""Settings for the Excel gene-name-error audit (Plan 1).

Everything that a reviewer might ask about lives here: years, seed, sample
targets, file-size caps, and the exact detection rules.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REF = DATA / "reference"        # gene lists, HGNC file, published rates
CAND = DATA / "candidates"      # candidate PMCID lists per year
DB_PATH = DATA / "scan.sqlite"  # record-level results (resumable)
OUT = ROOT / "output"           # tables, figures, summary

# ---- study design ---------------------------------------------------------
YEARS = list(range(2019, 2026))          # 2019 to 2025
SEED = 20261117                          # fixed seed, date the plan was written
TARGET_ELIGIBLE_PER_YEAR = 1500          # plan: about 1,500 eligible articles per year
BEFORE_YEARS = (2019, 2020)              # "before" period in the hypothesis
AFTER_YEARS = (2024, 2025)               # "after" period

# Same population as the 2021 study ("genom*" in title/abstract). NOTE: in the
# current PMC search, [Abstract] behaves like All Fields (139k hits for 2019
# instead of ~27k), so [Title/Abstract] is used. Checked 21 Nov 2026: this
# query contains 96.6% of the original 2019 article list and 96.5% of its
# articles with gene lists.
CANDIDATE_QUERY = "genom*[Title/Abstract] AND open_access[Filter]"

# ---- network --------------------------------------------------------------
NCBI_TOOL = "excel_gene_errors"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "hu.haoxi@northeastern.edu")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")   # optional, raises the rate limit
S3_BASE = "https://pmc-oa-opendata.s3.amazonaws.com"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/tsv/"
            "tsv/hgnc_complete_set.txt")
ORIG_RAW = "https://raw.githubusercontent.com/markziemann/GeneNameErrors2020/master"
ORIG_SPECIES_LISTS = [
    "Athaliana", "Celegans", "Dmelanogaster", "Drerio", "Ggallus", "Hsapiens",
    "Mmusculus", "Osativa", "Rnorvegicus", "Scerevisiae",
]
ORIG_YEARS = list(range(2014, 2021))
N_THREADS = int(os.environ.get("GNE_THREADS", "16"))

# ---- file handling --------------------------------------------------------
MAX_FILE_MB = 50            # bigger files are logged as "too_large", not silently dropped
MAX_ZIP_MEMBERS = 60
MAX_ZIP_UNCOMPRESSED_MB = 300
MAX_ROWS_PER_SHEET = 600_000
HEAD_ROWS = 100             # rows buffered to find gene columns
TOP_N_CELLS = 20            # original rule: look at the top 20 rows of a column
MIN_GENE_HITS = 5           # original rule: more than 4 of the top 20 are gene names
HEADER_ROWS_SKIPPED = 2     # original: header + first row deleted (sed 1,2d) before the regex

# Article types that are not research articles (plan: exclude reviews, corrections).
EXCLUDED_ARTICLE_TYPES = {
    "review-article", "systematic-review", "meta-analysis", "correction",
    "retraction", "editorial", "letter", "discussion", "news", "obituary",
    "expression-of-concern", "reply", "comment", "book-review",
    "article-commentary", "addendum", "other", "in-brief", "product-review",
}

# Journals named in the 2016 study (Ziemann, Eren, El-Osta, Genome Biology).
# All 18 journals from Table 1 of that paper (checked against the table on the
# Genome Biology site, 2026-11-23).
JOURNALS_2016 = [
    "Molecular Biology and Evolution", "Bioinformatics", "DNA Research",
    "Genome Biology and Evolution", "Nucleic Acids Research", "Genome Biology",
    "Nature Genetics", "Genome Research", "Genes & Development", "Nature",
    "Science", "PLoS ONE", "BMC Bioinformatics",
    "BMC Genomics", "Human Molecular Genetics", "RNA",
    "PLoS Computational Biology", "PLoS Biology",
]
HGNC_RENAMING_YEAR = 2020
EXCEL_UPDATE = 2023.8   # October 2023, as a decimal year for plotting
