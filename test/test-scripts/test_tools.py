"""
Installation test: verify each tool (a) runs without error on a small input and
(b) produces the SPECIFIC outputs this project depends on.

Writes a PASS/FAIL report to results/test-data/test_report.tsv.

Example (from bench/):
  python test/test-scripts/test_tools.py
  python test/test-scripts/test_tools.py -t fastoma orthophy (to test specific tools)
"""

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

# Scripts are in test/test-scripts/; config.py is in config/. Both are
# inside the project root (bench/), so config is two directories up.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config"))
import config

results = []

# Function to record one check result and print it
def record(tool, check, status, detail=""):
    results.append((tool, check, status, detail))
    print(f"[{status:4s}] {tool:<12} {check:<26} {detail}")

# Function to check whether a path has an output: a file must be non-empty,
# and a directory must contain at least one file. Several tools create their
# output folders before doing any work, so there just being a directory isn't enoguh.
def has_content(p):
    if os.path.isdir(p):
        return any(os.path.isfile(f) and os.path.getsize(f) > 0
                   for f in glob.glob(os.path.join(p, "**", "*"), recursive=True))
    return os.path.getsize(p) > 0

# Function to check that a pattern matches at least one path with content
def check_path(tool, label, pattern):
    hits = [p for p in glob.glob(pattern, recursive=True) if has_content(p)]
    if hits:
        record(tool, label, "PASS", os.path.relpath(hits[0]))
    else:
        record(tool, label, "FAIL", f"not found: {os.path.relpath(pattern)}")
    return bool(hits)

# Function to run a command in a tool's conda env, logging stdout and stderr
def run_in_env(tool, cmd, logfile, cwd=None):
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    with open(logfile, "w") as log:
        proc = subprocess.run(config.conda_run(tool) + cmd,
                              stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
    if proc.returncode == 0:
        record(tool, "runs", "PASS", "exit 0")
    else:
        record(tool, "runs", "FAIL",
               f"exit {proc.returncode}, see {os.path.basename(logfile)}")
    return proc.returncode == 0

# Function to copy every proteome file from one directory into another
def copy_proteomes(src, dest):
    for f in glob.glob(os.path.join(src, "*" + config.SHARED_EXT)):
        shutil.copy(f, dest)

# Function to write FastOMA's parameters as JSON
def write_fastoma_params(run_dir):
    params = {
        "input": "in_folder",
        "output_folder": "out",
        "omamer_db": "in_folder/omamerdb.h5",
        "species_tree": "in_folder/species_tree.nwk",
        "write_genetrees": True,    # gene trees are opt-in, off by default
        "write_msas": True,
    }
    with open(os.path.join(run_dir, "params.json"), "w") as fh:
        json.dump(params, fh, indent=2)

# Function to stage OrthoPhy input as a taxonomy-nested directory
# OrthoPhy reads taxonomy from directory nesting
def stage_orthophy_input(indir, staging):
    shutil.rmtree(staging, ignore_errors=True)
    tax = {}
    with open(config.TAXONOMY) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            sp, path = line.rstrip("\n").split("\t")
            tax[sp] = path
    staged, missing = 0, []
    for f in sorted(glob.glob(os.path.join(indir, "*" + config.SHARED_EXT))):
        sp = os.path.basename(f)[: -len(config.SHARED_EXT)]
        if sp not in tax:
            missing.append(sp)
            continue
        dest_dir = os.path.join(staging, tax[sp])
        os.makedirs(dest_dir, exist_ok=True)
        shutil.copy2(f, os.path.join(dest_dir, os.path.basename(f)))
        staged += 1
    if missing:
        raise ValueError(f"species missing from taxonomy.tsv: {missing}")
    return staged

# Function to test OrthoFinder
def test_orthofinder(indir, outdir):
    d = os.path.join(outdir, "orthofinder")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    copy_proteomes(indir, d)
    run_in_env("orthofinder",
               ["orthofinder", "-f", d, "-M", "msa", "-A", "mafft",
                "-T", config.OF_TREE_METHOD, "-t", str(config.THREADS)],
               os.path.join(config.LOG_DIR, "test_orthofinder.log"))
    res = os.path.join(d, "OrthoFinder", "Results_*")
    check_path("orthofinder", "orthogroups_table",
               os.path.join(res, "Orthogroups", "Orthogroups.tsv"))
    # v3 renamed Gene_Trees (v2) to Resolved_Gene_Trees: STRIDE-rooted trees
    check_path("orthofinder", "per_group_gene_trees",
               os.path.join(res, "Resolved_Gene_Trees"))
    check_path("orthofinder", "phylogenetic_hogs",
               os.path.join(res, "Phylogenetic_Hierarchical_Orthogroups"))
    check_path("orthofinder", "duplication_events",
               os.path.join(res, "Gene_Duplication_Events"))
    check_path("orthofinder", "sequence_id_map",
               os.path.join(res, "WorkingDirectory", "SequenceIDs.txt"))

# Function to test Broccoli
def test_broccoli(indir, outdir):
    d = os.path.join(outdir, "broccoli")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    run_in_env("broccoli",
               [config.BROCCOLI_CMD, "-dir", os.path.abspath(indir),
                "-o", os.path.abspath(d), "-threads", str(config.THREADS)],
               os.path.join(config.LOG_DIR, "test_broccoli.log"), cwd=d)
    check_path("broccoli", "orthologous_groups",
               os.path.join(d, "dir_step3", "table_OGs_protein_names.txt"))
    check_path("broccoli", "orthologous_pairs",
               os.path.join(d, "dir_step4", "orthologous_pairs.txt"))
    # Query-centric phylome, NOT one tree per orthogroup — different from OrthoFinder
    check_path("broccoli", "phylome_trees",
               os.path.join(d, "dir_step2", "phylomes.tsv"))

# Function to test FastOMA
def test_fastoma(indir, outdir):
    d = os.path.join(outdir, "fastoma")
    shutil.rmtree(d, ignore_errors=True)
    in_folder = os.path.join(d, "in_folder")
    os.makedirs(os.path.join(in_folder, "proteome"))
    copy_proteomes(indir, os.path.join(in_folder, "proteome"))
    for src, name, label in [(config.SPECIES_TREE, "species_tree.nwk", "species_tree_present"),
                             (config.OMAMER_DB, "omamerdb.h5", "omamer_db_present")]:
        if not os.path.isfile(src):
            record("fastoma", label, "FAIL", f"not found: {os.path.relpath(src)}")
            return
        # symlink, not copy: the database is 9.3 GB
        os.symlink(os.path.abspath(src), os.path.join(in_folder, name))
    write_fastoma_params(d)
    run_in_env("fastoma",
               ["nextflow", "run", "dessimozlab/FastOMA",
                "-r", config.FASTOMA_REVISION, "-params-file", "params.json"],
               os.path.join(config.LOG_DIR, "test_fastoma.log"), cwd=d)
    # RootHOGs = full families incl. paralogs. OrthologousGroups collapses to one
    # gene per species and would discard the paralogy under study.
    check_path("fastoma", "root_hogs", os.path.join(d, "out", "RootHOGs.tsv"))
    check_path("fastoma", "orthoxml", os.path.join(d, "out", "*.orthoxml"))
    check_path("fastoma", "gene_trees", os.path.join(d, "out", "**", "*enetree*"))
    check_path("fastoma", "msas", os.path.join(d, "out", "**", "msa*"))


# Function to test OrthoPhy
# --no_concatenate skips the concatenated species tree, which this project does
# not use and which fails on sparse test data. --force allows re-runs.
def test_orthophy(indir, outdir):
    d = os.path.join(outdir, "orthophy")
    shutil.rmtree(d, ignore_errors=True)
    if not os.path.isfile(config.ORTHOPHY_PY):
        record("orthophy", "runs", "SKIP", "not cloned into tools/OrthoPhy")
        return
    staging = os.path.join(outdir, "orthophy_input")
    try:
        n = stage_orthophy_input(indir, staging)
        record("orthophy", "taxonomy_staged", "PASS", f"{n} species nested by clade")
    except (ValueError, FileNotFoundError) as exc:
        record("orthophy", "taxonomy_staged", "FAIL", str(exc))
        return
    run_in_env("orthophy",
               ["python", config.ORTHOPHY_PY,
                "--indir", os.path.abspath(staging), "--out", os.path.abspath(d),
                "--threads", str(config.THREADS), "--no_concatenate", "--force"],
               os.path.join(config.LOG_DIR, "test_orthophy.log"))
    # Output layout is undocumented: OrthoPhy writes numbered step folders
    # (01_fasta, ..., 05_orthologs). List them with file counts so each run
    # records how far it got and where the gene trees land.
    if os.path.isdir(d):
        steps = sorted(e for e in os.listdir(d) if os.path.isdir(os.path.join(d, e)))
        summary = ", ".join(
            f"{e}({sum(1 for f in glob.glob(os.path.join(d, e, '**', '*'), recursive=True) if os.path.isfile(f))})"
            for e in steps)
        record("orthophy", "step_folders", "INFO", summary or "none")
    # Target the folder itself: a loose pattern like *rtholog* also matches the
    # tool's own log file (Ortholog-Finder2.log) and passes falsely.
    check_path("orthophy", "ortholog_sets", os.path.join(d, "05_orthologs"))
    # The check that decides whether OrthoPhy stays in the project
    check_path("orthophy", "gene_trees", os.path.join(d, "**", "*.tre*"))


TESTS = {
    "orthofinder": test_orthofinder,
    "broccoli": test_broccoli,
    "fastoma": test_fastoma,
    "orthophy": test_orthophy,
}


# Function to validate command line arguments
def check_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default=config.TEST_DIR,
                        help="directory of standardized test FASTAs "
                             "(default: test/test-data/prepped)")
    parser.add_argument("-o", "--output", default=config.TEST_RESULTS_DIR,
                        help="output directory (default: results/test-data)")
    parser.add_argument("-t", "--tools", nargs="+", default=config.CLI_TOOLS,
                        choices=config.CLI_TOOLS,
                        help="subset of tools to test (default: all)")
    args = parser.parse_args()
    if not os.path.isdir(args.input):
        parser.error(f"input directory not found: {args.input}")
    if not glob.glob(os.path.join(args.input, "*" + config.SHARED_EXT)):
        parser.error(f"no {config.SHARED_EXT} files in {args.input} — run prep_inputs.py first")
    return args


# Function to run the tests and write the report
def main(args):
    config.ensure_dirs()
    os.makedirs(args.output, exist_ok=True)
    for tool in args.tools:
        print(f"\n=== {tool} ===")
        TESTS[tool](args.input, args.output)
    report = os.path.join(args.output, "test_report.tsv")
    with open(report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tool", "check", "status", "detail"])
        w.writerows(results)
    fails = [r for r in results if r[2] == "FAIL"]
    print(f"\nReport: {os.path.relpath(report)}")
    print(f"{len(fails)} FAIL / {len(results)} checks")

if __name__ == "__main__":
    main(check_args())