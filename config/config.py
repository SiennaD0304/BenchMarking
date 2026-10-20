"""
Central configuration. Every script imports from here. No script hardcodes
paths or parameters.
"""

import os

# --- Project layout ---
# config.py stays in config/, one level below the project root (bench/).
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw_proteomes")     # downloaded, untouched
SHARED_DIR = os.path.join(DATA_DIR, "shared")         # standardized — all tools read this
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
LOG_DIR = os.path.join(RESULTS_DIR, "logs")
TOOLS_DIR = os.path.join(PROJECT_ROOT, "tools")       # source-distributed tools

# --- Installation test data ---
TEST_DATA_DIR = os.path.join(PROJECT_ROOT, "test", "test-data")
TEST_RAW_DIR = os.path.join(TEST_DATA_DIR, "test_raw")      # fetched, unprocessed
TEST_DIR = os.path.join(TEST_DATA_DIR, "prepped")           # standardized test set
TEST_RESULTS_DIR = os.path.join(RESULTS_DIR, "test-data")   # installation test output

# --- Shared input format ---
# Broccoli v2 reads only .fas/.fasta/.faa — NOT .fa. OrthoFinder and FastOMA
# accept all of these, so .faa is the extension every tool handles.
SHARED_EXT = ".faa"

# --- Tool parameters that must be pinned ---
# OrthoFinder's internal tree method (-T). iqtree3 fails on sparse alignments
# where a species contributes no single-copy orthogroups (all-gap rows), which
# happens on small test sets; fasttree tolerates them and matches the method
# Broccoli uses internally, removing one confound from the comparison.
OF_TREE_METHOD = "fasttree"      # fasttree | fasttree_fastest | raxml | raxml-ng | iqtree3

# FastOMA Nextflow revision. "main" is a moving development branch; pin a
# release tag (latest: v0.5.1) before benchmark runs, and check that it
# uses "input" rather than the older "input_folder".
FASTOMA_REVISION = "main"

# --- Run parameters ---
THREADS = 8               # fixed across tools so runtimes are comparable
REPLICATES = 3            # runs per condition, for variance
SCALE_SIZES = [5, 10, 15, 20, 25]   # species counts for the scaling curve
SEED = 12345              # pin stochastic tree search

# --- Auxiliary inputs (FastOMA needs a species tree and the OMAmer database,
#     OrthoPhy needs a taxonomy) ---
SPECIES_TREE = os.path.join(DATA_DIR, "species_tree.nwk")
TAXONOMY = os.path.join(DATA_DIR, "taxonomy.tsv")
# LUCA.h5 from omabrowser.org. Clade-restricted databases (Metazoa,
# Saccharomyceta) do not cover this taxon range. FastOMA expects it inside its
# input folder named omamerdb.h5; always pass --omamer_db explicitly, or
# FastOMA will download 9.9 GB on every run.
OMAMER_DB = os.path.join(DATA_DIR, "omamerdb.h5")

# --- Conda environment names (matches the `name:` field in environments/*.yml) ---
ENVS = {
    "orthofinder": "of3",
    "broccoli": "broccoli",
    "fastoma": "fastoma",
    "orthophy": "orthophy",
    "trees": "treepipe",
}

# --- Tools under test. OrthoVenn3 is web-based and handled manually. ---
CLI_TOOLS = ["orthofinder", "broccoli", "fastoma", "orthophy"]
ALL_TOOLS = CLI_TOOLS + ["orthovenn3"]

# --- Tool entry points ---
# Broccoli v2 is installed by bioconda.
BROCCOLI_CMD = "broccoli"
# OrthoPhy is installed from source (github.com/tmakwtnb/OrthoPhy).
ORTHOPHY_DIR = os.path.join(TOOLS_DIR, "OrthoPhy")
ORTHOPHY_PY = os.path.join(ORTHOPHY_DIR, "OrthoPhy.py")

# Function to create every project directory if it does not already exist
def ensure_dirs():
    for d in (RAW_DIR, SHARED_DIR, TEST_RAW_DIR, TEST_DIR,
              RESULTS_DIR, TEST_RESULTS_DIR, LOG_DIR, TOOLS_DIR):
        os.makedirs(d, exist_ok=True)

# Function to build a `conda run` command prefix for a given tool's environment
def conda_run(tool):
    return ["conda", "run", "-n", ENVS[tool], "--no-capture-output"]