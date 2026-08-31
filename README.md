# bpseudolongum assign spacers a unique identifier 
"""
assign_spacer_ids.py

Assign unique identifiers to CRISPR spacer sequences obtained from
CRISPRCasTyper output.

The input directory should contain one FASTA file per CRISPR array.
Each FASTA record represents one spacer, with headers in the format:

    >array_id:spacer_number

For example:

    >Cow_4_s37_scaffold_166416_1:1
    ACTGACCTGATCGATCGATCGATCGATCGATCGAT

    >Cow_4_s37_scaffold_166416_1:2
    TGCGATCGATCGTACGATCGATCGATCGATCGATC

The number following the colon indicates the position of the spacer
within the CRISPR array.

Each distinct spacer nucleotide sequence is assigned a global identifier
(s1, s2, s3, ...). Identical spacer sequences receive the same identifier,
including when the same spacer occurs in different CRISPR arrays or genomes.

Two tab-separated output files are generated:

    spacer_mapping.tsv
        Unique spacer ID and nucleotide sequence.

    array_index.tsv
        CRISPR array ID and the ordered spacer IDs in that array.

Usage:

    python assign_spacer_ids.py input_directory/ -o output_directory/

Requires Python 3 and uses only modules from the Python standard library.
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict


def read_fasta(path):
    """Read a FASTA file and yield header, sequence pairs."""

    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()

                header = line[1:].strip()
                sequence = []

            else:
                sequence.append(line)

        if header is not None:
            yield header, "".join(sequence).upper()


def parse_header(header):
    """Return the CRISPR array ID and spacer position from a FASTA header."""

    match = re.match(r"^(.*):(\d+)$", header)

    if match is None:
        raise ValueError(
            f"Unexpected FASTA header: {header}\n"
            "Expected format: >array_id:spacer_number"
        )

    array_id = match.group(1)
    spacer_number = int(match.group(2))

    return array_id, spacer_number


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Assign unique IDs to CRISPR spacer sequences and record "
            "their order within each CRISPR array."
        )
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing CRISPRCasTyper spacer FASTA files."
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Directory where output files will be written."
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(
        path for path in input_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() in {".fa", ".fasta", ".fna"}
    )

    if not fasta_files:
        raise SystemExit(f"No FASTA files found in {input_dir}")

    arrays = defaultdict(list)

    for fasta_file in fasta_files:
        for header, sequence in read_fasta(fasta_file):

            array_id, spacer_number = parse_header(header)

            arrays[array_id].append(
                (spacer_number, sequence)
            )

    sequence_to_id = {}
    id_to_sequence = []

    for array_id in sorted(arrays):

        ordered_spacers = sorted(
            arrays[array_id],
            key=lambda x: x[0]
        )

        for _, sequence in ordered_spacers:

            if sequence not in sequence_to_id:

                spacer_id = f"s{len(sequence_to_id) + 1}"

                sequence_to_id[sequence] = spacer_id
                id_to_sequence.append(
                    (spacer_id, sequence)
                )

    mapping_file = outdir / "spacer_mapping.tsv"

    with open(mapping_file, "w") as handle:

        handle.write("spacer_id\tsequence\n")

        for spacer_id, sequence in id_to_sequence:
            handle.write(
                f"{spacer_id}\t{sequence}\n"
            )

    array_file = outdir / "array_index.tsv"

    with open(array_file, "w") as handle:

        handle.write("array_id\tspacer_ids\n")

        for array_id in sorted(arrays):

            ordered_spacers = sorted(
                arrays[array_id],
                key=lambda x: x[0]
            )

            spacer_ids = [
                sequence_to_id[sequence]
                for _, sequence in ordered_spacers
            ]

            handle.write(
                f"{array_id}\t{','.join(spacer_ids)}\n"
            )

    print(f"FASTA files processed: {len(fasta_files)}")
    print(f"CRISPR arrays identified: {len(arrays)}")
    print(f"Unique spacer sequences: {len(sequence_to_id)}")
    print(f"Spacer mapping: {mapping_file}")
    print(f"Array index: {array_file}")


if __name__ == "__main__":
    main()
