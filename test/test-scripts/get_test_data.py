"""
Pulls the entire installation test dataset from NCBI, so nothing needs to be
downloaded manually.

WHY THE INSTALLATION TEST ISN'T SIMPLER
A single gene with one copy per species does not ensure that every part of every
tool works. With no duplication nodes, Broccoli may build no trees at all,
OrthoFinder's STRIDE rooting has nothing to root on, and FastOMA's RootHOGs.tsv
becomes identical to OrthologousGroups.tsv — so a key distinction for this project
would be untested. The set therefore needs paralogy and more than one family.

WHAT IS FETCHED FROM NCBI
  1. Seven focal accessions, by accession number.
  2. Four gene families, by gene symbol and organism, as a difficulty gradient:
     phosducin    Focal family. Includes PDC/PDCL/PDCL2/PDCL3/TXNDC9 plus the
                  yeast paralogs PLP1/PLP2. Its ancient duplications and quickly
                  evolving Plasmodium member make it the hardest case.
     thioredoxin  Difficult decoy. TXN/TXN2/TXNL1 have the same thioredoxin
                  domain as TXNDC9, so a tool that groups them with the
                  phosducin family is making a significant error.
     14-3-3       Simple decoy. The YWHA family occurs several times in every
                  species in this set and is unrelated to the focal family.
     actin        Simple decoy. Also present several times in every species,
                  but so conserved that tools should cluster it easily.

Searching by symbol rather than by a hardcoded accession list means results
follow current RefSeq data rather than going out of date.

Every header is annotated with [family=...] so that if a tool misgroups
something, the report shows which family it pulled in.

Requires network access and an email address (NCBI E-utilities policy).

EXAMPLE USAGE (from bench/)
  python test/test-scripts/get_test_data.py -e user@email.com
"""

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request

# Scripts live in test/test-scripts/; config.py lives in config/. Both sit
# inside the project root (bench/), so config is two directories up.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config"))
import config

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Species tag -> NCBI organism name. Tags become the species prefix in each ID.
SPECIES = {
    "Hsap": "Homo sapiens",
    "Drer": "Danio rerio",
    "Asag": "Anolis sagrei",
    "Dmel": "Drosophila melanogaster",
    "Cele": "Caenorhabditis elegans",
    "Scer": "Saccharomyces cerevisiae",
    "Pber": "Plasmodium berghei",
}

# The seven focal sequences, fetched by accession so they are always present
# even where symbol search returns nothing (Anolis and Plasmodium are thinly
# annotated and may not respond to [Gene Name] queries).
SEED_ACCESSIONS = {
    "Hsap": ["NP_005774.2"],       # TXNDC9
    "Drer": ["NP_001315504.1"],    # TXNDC9
    "Asag": ["XP_060625760.1"],    # TXNDC9 (truncated record)
    "Dmel": ["Q9VGV8.1"],          # Phosducin-like protein 3
    "Cele": ["NP_498410.1"],       # TXNDC9
    "Scer": ["NP_010469.3"],       # Plp1p
    "Pber": ["XP_034422726.1"],    # Phosducin-like protein 3
}

# Family -> gene symbols to search. Symbols differ between lineages, so each
# family lists aliases; absent genes simply return nothing.
FAMILIES = {
    "phosducin": ["PDC", "PDCL", "PDCL2", "PDCL3", "TXNDC9",
                  "PLP1", "PLP2", "phosducin-like"],
    "thioredoxin": ["TXN", "TXN2", "TXNL1", "TRX1", "TRX2"],
    "14-3-3": ["YWHAB", "YWHAE", "YWHAZ", "YWHAG", "YWHAH", "YWHAQ",
               "SFN", "BMH1", "BMH2", "14-3-3"],
    "actin": ["ACTB", "ACTG1", "ACT1", "ACT5C", "act-1"],
}

# Function to call NCBI E-utilities with retries and rate limiting
def eutils(endpoint, params, email, retries=3):
    params = dict(params)
    params["email"] = email
    params["tool"] = "bioi500-ortho-benchmark"
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return fh.read().decode()
        except Exception as exc:
            if attempt == retries - 1:
                sys.stderr.write(f"[warn] {endpoint} failed: {exc}\n")
                return ""
            time.sleep(2 ** attempt)
        finally:
            time.sleep(0.4)      # NCBI: max 3 requests/sec without an API key
    return ""

# Function to search the protein database for one gene symbol in one organism
def search_protein(symbol, organism, email, retmax):
    # refseq filter keeps results curated; without it a search returns many
    # near-duplicate submissions that inflate the test set
    query = f'{symbol}[Gene Name] AND "{organism}"[Organism] AND refseq[filter]'
    xml = eutils("esearch", {"db": "protein", "term": query,
                             "retmax": retmax, "retmode": "xml"}, email)
    return [line.strip()[4:-5] for line in xml.splitlines()
            if line.strip().startswith("<Id>")]

# Function to download protein sequences in FASTA format by UID or accession
def fetch_fasta(ids, email):
    if not ids:
        return ""
    return eutils("efetch", {"db": "protein", "id": ",".join(ids),
                             "rettype": "fasta", "retmode": "text"}, email)

# Function to split FASTA into (accession, record) pairs
def split_records(fasta):
    out = []
    for block in fasta.split("\n>"):
        block = block.strip()
        if not block:
            continue
        if not block.startswith(">"):
            block = ">" + block
        out.append((block[1:].split()[0], block))
    return out

# Function to tag a record's header with the family it came from
def tag_family(record, family):
    header, _, body = record.partition("\n")
    return f"{header} [family={family}]\n{body}"

# Function to collect every record for one species
def fetch_species(tag, organism, families, email, max_per_gene):
    records, seen = [], set()
    for acc, record in split_records(fetch_fasta(SEED_ACCESSIONS[tag], email)):
        if acc not in seen:
            seen.add(acc)
            records.append(tag_family(record, "phosducin"))
    seed_n = len(records)
    counts = {}
    for family in families:
        before = len(records)
        for symbol in FAMILIES[family]:
            uids = search_protein(symbol, organism, email, max_per_gene)
            if not uids:
                continue
            for acc, record in split_records(fetch_fasta(uids, email)):
                if acc in seen:
                    continue
                seen.add(acc)
                records.append(tag_family(record, family))
        counts[family] = len(records) - before
    return records, seed_n, counts

# Function to retrieve command line arguments
def check_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", default=config.TEST_RAW_DIR,
                        help="directory for per-species FASTA "
                             "(default: test/test-data/test_raw)")
    parser.add_argument("-e", "--email", required=True,
                        help="your email — required by NCBI E-utilities policy")
    parser.add_argument("-f", "--families", nargs="+", default=list(FAMILIES),
                        choices=list(FAMILIES),
                        help="subset of families to fetch (default: all)")
    parser.add_argument("--max-per-gene", type=int, default=5,
                        help="cap hits per gene symbol per species (default: 5)")
    return parser.parse_args()

# Function to fetch every species and write per-species FASTA files
def main(args):
    config.ensure_dirs()
    os.makedirs(args.output, exist_ok=True)
    header = f"{'species':<8} {'seed':>5}" + "".join(
        f" {f:>12}" for f in args.families) + f" {'total':>7}"
    print(header)
    print("-" * len(header))
    total = 0
    for tag, organism in SPECIES.items():
        records, seed_n, counts = fetch_species(
            tag, organism, args.families, args.email, args.max_per_gene)
        with open(os.path.join(args.output, f"{tag}.fa"), "w") as fh:
            fh.write("\n".join(records) + "\n")
        row = f"{tag:<8} {seed_n:>5}" + "".join(
            f" {counts[f]:>12}" for f in args.families) + f" {len(records):>7}"
        print(row)
        total += len(records)
    print("-" * len(header))
    print(f"{total} sequences written to {args.output}")
    if total < 40:
        print("\n[warn] Fewer sequences than expected. Anolis sagrei and "
              "Plasmodium berghei are thinly annotated and may not respond to "
              "gene-symbol search; consider substituting Anolis carolinensis, "
              "or adding accessions to SEED_ACCESSIONS by hand.")

if __name__ == "__main__":
    main(check_args())