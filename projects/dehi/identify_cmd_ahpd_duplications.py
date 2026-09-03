"""Identify tandemly duplicated CMD-core architectures with a local profile HMM.

The seed profile is built from the C-terminal CMD-core region of compact,
explicitly CMD-annotated PF02627 proteins. Longer CMD-like and AhpD proteins are retained as duplication
outcomes only when the seed HMM detects two significant, substantially covered,
non-overlapping matches. The intact double-core proteins are used for later
comparison with the intact DehI halves; extracted core matches are written to
validate the architecture and boundary.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyhmmer


ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
CMD_INPUT = RESULTS / "uniprot_PF02627_nr60.fasta"
AHPD_INPUT = RESULTS / "uniprot_ahpd_IPR004674_nr60.fasta"
SEED_FASTA = RESULTS / "cmd_single_core_seed.fasta"
SINGLE_MONOMER_FASTA = RESULTS / "cmd_single_monomer_80_119.fasta"
SEED_MSA = RESULTS / "cmd_single_core_seed_msa.fasta"
SEED_A3M = RESULTS / "cmd_single_core_seed.a3m"
SEED_STOCKHOLM = RESULTS / "cmd_single_core_seed.sto"
SEED_HMM = RESULTS / "cmd_single_core_seed.hmm"
CMD_CANDIDATES = RESULTS / "cmd_double_core_candidates.fasta"
AHPD_CANDIDATES = RESULTS / "ahpd_double_core_candidates.fasta"
ARCHITECTURE_REPORT = RESULTS / "cmd_ahpd_duplication_architectures.tsv"
NONSTANDARD = str.maketrans({character: "X" for character in "BJOUZ*"})

# The seed represents compact, single-copy-like PF02627 proteins. Candidate
# bounds are intentionally broader; architecture is decided by two HMM hits.
SEED_MIN_LENGTH = 80
SEED_MAX_LENGTH = 119
CANDIDATE_MIN_LENGTH = 120
CANDIDATE_MAX_LENGTH = 210
COMPACT_DUPLICATION_MAX_LENGTH = 150
MAX_DOMAIN_I_EVALUE = 1e-3
MIN_HMM_COVERAGE = 0.45
MAX_DOMAIN_OVERLAP = 10
SEED_CORE_START_FRACTION = 0.45


@dataclass(frozen=True)
class CoreHit:
    start: int
    end: int
    hmm_start: int
    hmm_end: int
    score: float
    i_evalue: float

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def normalized_copy(record: SeqRecord) -> SeqRecord:
    copied = record[:]
    copied.seq = Seq(str(record.seq).upper().translate(NONSTANDARD))
    return copied


def is_cmd_annotation(record: SeqRecord) -> bool:
    description = record.description.lower()
    return (
        "carboxymuconolactone decarboxylase" in description
        and "ahpd" not in description
    )


def prepare_inputs() -> None:
    cmd = [normalized_copy(record) for record in SeqIO.parse(CMD_INPUT, "fasta")]
    ahpd = [normalized_copy(record) for record in SeqIO.parse(AHPD_INPUT, "fasta")]

    compact_seed_proteins = [
        record
        for record in cmd
        if is_cmd_annotation(record)
        and SEED_MIN_LENGTH <= len(record.seq) <= SEED_MAX_LENGTH
    ]
    # Structural comparisons define the CMD core as helices alpha4-alpha6,
    # located in the C-terminal portion of compact six-helix CMD-family
    # proteins. Retaining the final 55% avoids training on the variable
    # N-terminal helices while allowing proportional length variation.
    seed = []
    for record in compact_seed_proteins:
        start = round(len(record.seq) * SEED_CORE_START_FRACTION)
        seed.append(
            SeqRecord(
                record.seq[start:],
                id=f"{record.id}|core_seed|{start + 1}-{len(record.seq)}",
                description=record.description,
            )
        )
    cmd_candidates = [
        record
        for record in cmd
        if is_cmd_annotation(record)
        and CANDIDATE_MIN_LENGTH <= len(record.seq) <= CANDIDATE_MAX_LENGTH
    ]
    ahpd_candidates = [
        record
        for record in ahpd
        if CANDIDATE_MIN_LENGTH <= len(record.seq) <= CANDIDATE_MAX_LENGTH
    ]

    SeqIO.write(seed, SEED_FASTA, "fasta")
    SeqIO.write(compact_seed_proteins, SINGLE_MONOMER_FASTA, "fasta")
    SeqIO.write(cmd_candidates, CMD_CANDIDATES, "fasta")
    SeqIO.write(ahpd_candidates, AHPD_CANDIDATES, "fasta")
    write_midpoint_cohort(
        "cmd_compact_120_150",
        [record for record in cmd_candidates if len(record.seq) <= COMPACT_DUPLICATION_MAX_LENGTH],
    )
    write_midpoint_cohort(
        "ahpd_compact_120_150",
        [record for record in ahpd_candidates if len(record.seq) <= COMPACT_DUPLICATION_MAX_LENGTH],
    )
    print(f"Single-core seed: {len(seed)} structure-guided fragments")
    print(f"CMD-like candidates: {len(cmd_candidates)} sequences")
    print(f"AhpD candidates: {len(ahpd_candidates)} sequences")


def write_midpoint_cohort(stem: str, records: list[SeqRecord]) -> None:
    """Write intact and midpoint-split sequences for a length-defined cohort."""
    first_repeats = []
    second_repeats = []
    for record in records:
        boundary = len(record.seq) // 2
        first_repeats.append(
            SeqRecord(
                record.seq[:boundary],
                id=f"{record.id}|midpoint_repeat1|1-{boundary}",
                description=record.description,
            )
        )
        second_repeats.append(
            SeqRecord(
                record.seq[boundary:],
                id=f"{record.id}|midpoint_repeat2|{boundary + 1}-{len(record.seq)}",
                description=record.description,
            )
        )
    SeqIO.write(records, RESULTS / f"{stem}.fasta", "fasta")
    SeqIO.write(first_repeats, RESULTS / f"{stem}_repeat1.fasta", "fasta")
    SeqIO.write(second_repeats, RESULTS / f"{stem}_repeat2.fasta", "fasta")
    print(f"{stem}: {len(records)} sequences")


def build_seed_hmm():
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(
        SEED_STOCKHOLM, digital=True, alphabet=alphabet
    ) as msa_file:
        msa = msa_file.read()
    if msa is None:
        raise ValueError(f"No alignment found in {SEED_MSA}")
    msa.name = b"CMD_single_core_seed"
    # reformat.pl marks the <50%-gap columns in the Stockholm RF line. Hand
    # architecture makes HMMER use exactly those states instead of treating
    # sparse insertion columns as additional core positions.
    builder = pyhmmer.plan7.Builder(alphabet, architecture="hand")
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    with SEED_HMM.open("wb") as handle:
        hmm.write(handle)
    return alphabet, hmm


def choose_two_domains(hit, hmm_length: int) -> tuple[CoreHit, CoreHit] | None:
    domains = []
    for domain in hit.domains:
        alignment = domain.alignment
        hmm_coverage = (
            alignment.hmm_to - alignment.hmm_from + 1
        ) / hmm_length
        if domain.i_evalue <= MAX_DOMAIN_I_EVALUE and hmm_coverage >= MIN_HMM_COVERAGE:
            domains.append(
                CoreHit(
                    start=alignment.target_from,
                    end=alignment.target_to,
                    hmm_start=alignment.hmm_from,
                    hmm_end=alignment.hmm_to,
                    score=domain.score,
                    i_evalue=domain.i_evalue,
                )
            )

    best = None
    best_score = float("-inf")
    for index, first in enumerate(domains):
        for second in domains[index + 1 :]:
            left, right = sorted((first, second), key=lambda domain: domain.start)
            overlap = left.end - right.start + 1
            if overlap <= MAX_DOMAIN_OVERLAP:
                score = left.score + right.score
                if score > best_score:
                    best = (left, right)
                    best_score = score
    return best


def scan_group(group: str, path: Path, alphabet, hmm):
    records = list(SeqIO.parse(path, "fasta"))
    by_id = {record.id: record for record in records}
    with pyhmmer.easel.SequenceFile(
        path, digital=True, alphabet=alphabet
    ) as sequence_file:
        sequences = list(sequence_file)

    (top_hits,) = pyhmmer.hmmer.hmmsearch(
        hmm,
        sequences,
        cpus=0,
        E=100.0,
        domE=100.0,
        incE=100.0,
        incdomE=100.0,
    )
    hits_by_name = {str(hit.name): hit for hit in top_hits}
    selected = []
    rows = []
    core_1 = []
    core_2 = []

    for record in records:
        hit = hits_by_name.get(record.id)
        pair = choose_two_domains(hit, hmm.M) if hit is not None else None
        if pair is None:
            rows.append((group, record.id, len(record.seq), 0, "", "", "", ""))
            continue

        first, second = pair
        selected.append(record)
        core_1.append(
            SeqRecord(
                record.seq[first.start - 1 : first.end],
                id=f"{record.id}|core1|{first.start}-{first.end}",
                description=record.description,
            )
        )
        core_2.append(
            SeqRecord(
                record.seq[second.start - 1 : second.end],
                id=f"{record.id}|core2|{second.start}-{second.end}",
                description=record.description,
            )
        )
        rows.append(
            (
                group,
                record.id,
                len(record.seq),
                1,
                f"{first.start}-{first.end}",
                f"{first.i_evalue:.3g}",
                f"{second.start}-{second.end}",
                f"{second.i_evalue:.3g}",
            )
        )

    stem = "cmd" if group == "CMD-like" else "ahpd"
    SeqIO.write(selected, RESULTS / f"{stem}_double_core.fasta", "fasta")
    SeqIO.write(core_1, RESULTS / f"{stem}_double_core_repeat1.fasta", "fasta")
    SeqIO.write(core_2, RESULTS / f"{stem}_double_core_repeat2.fasta", "fasta")
    print(f"{group}: selected {len(selected)} of {len(records)} candidates")
    return rows


def scan_candidates() -> None:
    alphabet, hmm = build_seed_hmm()
    rows = []
    rows.extend(scan_group("CMD-like", CMD_CANDIDATES, alphabet, hmm))
    rows.extend(scan_group("AhpD", AHPD_CANDIDATES, alphabet, hmm))
    with ARCHITECTURE_REPORT.open("w") as handle:
        handle.write(
            "group\taccession\tprotein_length\tdouble_core\t"
            "repeat1\trepeat1_i_evalue\trepeat2\trepeat2_i_evalue\n"
        )
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Architecture report: {ARCHITECTURE_REPORT}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage", choices=("prepare", "scan"), help="prepare FASTA files or scan candidates"
    )
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare_inputs()
    else:
        scan_candidates()


if __name__ == "__main__":
    main()
