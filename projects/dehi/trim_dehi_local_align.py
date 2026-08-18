from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


INPUT_FASTA = Path("results/dehi_msa.fasta")  # original unaligned sequences
OUTPUT_FASTA = Path("results/dehi_local_core_trimmed.fasta")
REPORT_FILE = Path("results/dehi_local_core_trim_report.tsv")

REFERENCE_TEXT = "3BJX"

# Require the local match to cover this fraction of 3BJX.
MIN_REFERENCE_COVERAGE = 0.70


records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

reference_record = next(
    (
        record
        for record in records
        if REFERENCE_TEXT.lower() in record.id.lower()
        or REFERENCE_TEXT.lower() in record.description.lower()
    ),
    None,
)

if reference_record is None:
    raise ValueError(f"Could not find a sequence containing {REFERENCE_TEXT!r}")

reference_sequence = str(reference_record.seq).replace("-", "").upper()


aligner = PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

# Reasonable protein-alignment gap penalties
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


trimmed_records = []
report_rows = []

for record in records:
    sequence = str(record.seq).replace("-", "").upper()

    if record.id == reference_record.id:
        trimmed_records.append(
            SeqRecord(
                Seq(sequence),
                id=record.id,
                description="reference_3BJX",
            )
        )
        report_rows.append(
            (
                record.id,
                len(sequence),
                1,
                len(sequence),
                len(sequence),
                1.0,
                "reference",
            )
        )
        continue

    alignments = aligner.align(reference_sequence, sequence)

    if len(alignments) == 0:
        report_rows.append(
            (record.id, len(sequence), "", "", "", 0.0, "no_alignment")
        )
        continue

    best_alignment = alignments[0]

    # aligned[0] contains blocks in the 3BJX reference.
    # aligned[1] contains corresponding blocks in the query protein.
    reference_blocks = best_alignment.aligned[0]
    query_blocks = best_alignment.aligned[1]

    reference_start = int(reference_blocks[0][0])
    reference_end = int(reference_blocks[-1][1])

    query_start = int(query_blocks[0][0])
    query_end = int(query_blocks[-1][1])

    reference_coverage = (
        reference_end - reference_start
    ) / len(reference_sequence)

    if reference_coverage < MIN_REFERENCE_COVERAGE:
        report_rows.append(
            (
                record.id,
                len(sequence),
                query_start + 1,
                query_end,
                query_end - query_start,
                round(reference_coverage, 3),
                "rejected_low_coverage",
            )
        )
        continue

    trimmed_sequence = sequence[query_start:query_end]

    trimmed_records.append(
        SeqRecord(
            Seq(trimmed_sequence),
            id=record.id,
            description=(
                f"local_3BJX_core "
                f"original_length={len(sequence)} "
                f"original_positions={query_start + 1}-{query_end} "
                f"reference_coverage={reference_coverage:.3f}"
            ),
        )
    )

    report_rows.append(
        (
            record.id,
            len(sequence),
            query_start + 1,
            query_end,
            len(trimmed_sequence),
            round(reference_coverage, 3),
            "retained",
        )
    )


SeqIO.write(trimmed_records, OUTPUT_FASTA, "fasta")

with REPORT_FILE.open("w") as handle:
    handle.write(
        "sequence_id\toriginal_length\tcore_start\tcore_end\t"
        "trimmed_length\treference_coverage\tstatus\n"
    )

    for row in report_rows:
        handle.write("\t".join(map(str, row)) + "\n")


print(f"Reference: {reference_record.id}")
print(f"Reference length: {len(reference_sequence)} aa")
print(f"Retained sequences: {len(trimmed_records)}")
print(f"Trimmed FASTA: {OUTPUT_FASTA}")
print(f"Report: {REPORT_FILE}")