"""
Standardize raw proteomes into the shared input directory.

Every tool reads from the same shared directory so "identical input" is
guaranteed.

Transformations applied to each proteome:
  1. truncate headers at first whitespace, strip database prefixes
     (sp|ACC|NAME -> ACC), and prefix with the species tag
  2. strip stop codons (*) — these break DIAMOND in some tools
  3. replace non-standard residues (U, O) with X
  4. drop duplicate IDs within a species
  5. drop sequences shorter than --min-length

Output files use the .faa extension: Broccoli v2 reads only .fas/.fasta/.faa,
not .fa, and .faa is accepted by every tool in this project.

Run in any environment with Python 3 (no external dependencies).

Example (from bench/):
  python test/test-scripts/prep_inputs.py -i test/test-data/test_raw -o test/test-data/prepped
"""

import argparse
import gzip
import os
import sys

# Scripts are in test/test-scripts/; config.py is in config/. Both are
# inside project root (bench/), so config is two directories up.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config"))
import config

FASTA_EXT = (".fa", ".faa", ".fasta", ".fas", ".pep")

# Function to open a plain or gzipped file
def smart_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

# Function to derive a species tag from a filename
def species_from_filename(path):
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for ext in FASTA_EXT:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name

# Function to reduce a sequence ID to a single bare accession
# SwissProt headers (sp|Q9VGV8.1|PHLP3_DROME) contain no whitespace, so
# truncating at the first space leaves embedded pipes. The standard format is
# Species|accession, so the extra fields must be removed or the ID ends up with four
# pipe-delimited parts instead of two.
def normalize_id(acc):
    if acc.startswith(("sp|", "tr|", "gi|", "ref|", "pdb|")):
        parts = acc.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return acc.replace("|", "_")

# Function to read a FASTA file, returning (header_id, sequence) pairs
def read_fasta(path):
    header, chunks = None, []
    with smart_open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = normalize_id(line[1:].split()[0])   # truncate, then strip db prefixes
                chunks = []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

# Function to clean a sequence string
def clean_sequence(seq):
    seq = seq.replace("*", "")                  # remove trailing and internal stops
    return seq.upper().replace("U", "X").replace("O", "X")

# Function to standardize a proteome and write it to the shared directory
def prep_one(infile, outdir, min_length):
    species = species_from_filename(infile)
    outfile = os.path.join(outdir, species + config.SHARED_EXT)
    seen, written, skipped_dup, skipped_short = set(), 0, 0, 0
    with open(outfile, "w") as out:
        for acc, seq in read_fasta(infile):
            if acc in seen:
                skipped_dup += 1
                continue
            seen.add(acc)
            seq = clean_sequence(seq)
            if len(seq) < min_length:
                skipped_short += 1
                continue
            # Standard ID format: Species|accession
            out.write(f">{species}|{acc}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            written += 1
    return {
        "species": species,
        "outfile": outfile,
        "written": written,
        "dup": skipped_dup,
        "short": skipped_short,
    }

# Function to validate command line arguments
def check_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default=config.RAW_DIR,
                        help="directory of raw proteome FASTAs (default: data/raw_proteomes)")
    parser.add_argument("-o", "--output", default=config.SHARED_DIR,
                        help="shared output directory (default: data/shared)")
    parser.add_argument("--min-length", type=int, default=30,
                        help="drop sequences shorter than this (default: 30)")
    args = parser.parse_args()
    if not os.path.isdir(args.input):
        parser.error(f"input directory not found: {args.input}")
    return args

# Function to standardize each proteome
def main(args):
    os.makedirs(args.output, exist_ok=True)
    infiles = sorted(
        os.path.join(args.input, f) for f in os.listdir(args.input)
        if f.endswith(FASTA_EXT) or f.endswith(tuple(e + ".gz" for e in FASTA_EXT))
    )
    if not infiles:
        sys.exit(f"No FASTA files found in {args.input}")
    print(f"{'species':<20} {'kept':>8} {'dup':>6} {'short':>6}")
    print("-" * 44)
    for infile in infiles:
        r = prep_one(infile, args.output, args.min_length)
        print(f"{r['species']:<20} {r['written']:>8} {r['dup']:>6} {r['short']:>6}")
    print("-" * 44)
    print(f"{len(infiles)} species written to {args.output}")

if __name__ == "__main__":
    main(check_args())