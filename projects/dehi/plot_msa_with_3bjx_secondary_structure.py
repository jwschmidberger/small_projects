"""Plot a wrapped MSA with the experimental 3BJX secondary structure on top.

The default reads deposited HELIX/SHEET records and needs only uv-installable Python
packages. Pass ``--source dssp`` to recalculate assignments with ``mkdssp``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlretrieve

from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from matplotlib.patches import FancyArrow, FancyBboxPatch
from pymsaviz import MsaViz


HELIX_CODES = frozenset("HGI")
STRAND_CODES = frozenset("EB")


def download_3bjx(destination: Path) -> Path:
    """Download 3BJX coordinates when they are not already available."""
    if not destination.exists():
        print(f"Downloading 3BJX coordinates to {destination}")
        urlretrieve("https://files.rcsb.org/download/3BJX.pdb", destination)
    return destination


def dssp_sequence_and_profile(pdb_file: Path, chain_id: str) -> tuple[str, str]:
    """Return the coordinate-derived sequence and DSSP profile for one chain."""
    from Bio.PDB.DSSP import DSSP

    structure = PDBParser(QUIET=True).get_structure("3BJX", pdb_file)
    model = structure[0]
    dssp = DSSP(model, str(pdb_file), dssp="mkdssp")

    amino_acids: list[str] = []
    assignments: list[str] = []
    for key in dssp.keys():
        if key[0] != chain_id:
            continue
        amino_acid = dssp[key][1].upper()
        amino_acids.append(amino_acid if amino_acid.isalpha() else "X")
        assignments.append(dssp[key][2] or "-")

    if not amino_acids:
        available = sorted({key[0] for key in dssp.keys()})
        raise ValueError(
            f"Chain {chain_id!r} has no DSSP residues; available chains: {available}"
        )
    return "".join(amino_acids), "".join(assignments)


def deposited_sequence_and_profile(pdb_file: Path, chain_id: str) -> tuple[str, str]:
    """Read wwPDB HELIX/SHEET assignments without an external executable."""
    structure = PDBParser(QUIET=True).get_structure("3BJX", pdb_file)
    chain = structure[0][chain_id]
    residues = [residue for residue in chain if is_aa(residue, standard=False)]
    sequence = "".join(
        seq1(residue.get_resname(), custom_map={"MSE": "M"}) for residue in residues
    )
    profile = ["-"] * len(residues)

    # Values are (kind, start residue number, end residue number). 3BJX has no
    # insertion-code ambiguity in its deposited secondary-structure ranges.
    ranges: list[tuple[str, int, int]] = []
    with pdb_file.open() as handle:
        for line in handle:
            if line.startswith("HELIX "):
                if line[19].strip() == chain_id and line[31].strip() == chain_id:
                    ranges.append(("H", int(line[21:25]), int(line[33:37])))
            elif line.startswith("SHEET "):
                if line[21].strip() == chain_id and line[32].strip() == chain_id:
                    ranges.append(("E", int(line[22:26]), int(line[33:37])))

    if not ranges:
        raise ValueError(f"No deposited HELIX/SHEET records found for chain {chain_id!r}")
    for kind, start_number, end_number in ranges:
        for index, residue in enumerate(residues):
            residue_number = residue.id[1]
            if start_number <= residue_number <= end_number:
                profile[index] = kind
    return sequence, "".join(profile)


def find_reference(alignment, keyword: str):
    # Prefer an exact FASTA record ID. Other records may mention 3BJX in their
    # descriptions (for example ``local_3BJX_core``) without being the reference.
    exact_matches = [record for record in alignment if record.id.lower() == keyword.lower()]
    if len(exact_matches) == 1:
        return exact_matches[0]

    matches = [
        record
        for record in alignment
        if record.description.lower().startswith(keyword.lower())
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one MSA record containing {keyword!r}; found {len(matches)}"
        )
    return matches[0]


def map_dssp_to_alignment(
    aligned_reference: str, dssp_sequence: str, dssp_profile: str
) -> tuple[str, float]:
    """Map DSSP assignments onto MSA columns through a sequence alignment."""
    reference_sequence = aligned_reference.replace("-", "").replace(".", "")
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    mapping_alignment = aligner.align(reference_sequence, dssp_sequence)[0]

    dssp_to_reference: dict[int, int] = {}
    matches = 0
    paired = 0
    for (ref_start, ref_end), (dssp_start, dssp_end) in zip(
        mapping_alignment.aligned[0], mapping_alignment.aligned[1]
    ):
        block_length = min(ref_end - ref_start, dssp_end - dssp_start)
        for offset in range(block_length):
            ref_index = ref_start + offset
            dssp_index = dssp_start + offset
            dssp_to_reference[dssp_index] = ref_index
            paired += 1
            matches += reference_sequence[ref_index] == dssp_sequence[dssp_index]

    reference_profile = ["-"] * len(reference_sequence)
    for dssp_index, ref_index in dssp_to_reference.items():
        reference_profile[ref_index] = dssp_profile[dssp_index]

    msa_profile: list[str] = []
    residue_index = 0
    for character in aligned_reference:
        if character in "-.":
            msa_profile.append("-")
        else:
            msa_profile.append(reference_profile[residue_index])
            residue_index += 1

    identity = matches / paired if paired else 0.0
    return "".join(msa_profile), identity


def runs(profile: str, accepted_codes: frozenset[str], start: int, end: int):
    """Yield half-open contiguous runs within one wrapped alignment panel."""
    run_start = None
    for column in range(start, end + 1):
        accepted = column < end and profile[column] in accepted_codes
        if accepted and run_start is None:
            run_start = column
        elif not accepted and run_start is not None:
            yield run_start, column
            run_start = None


def add_secondary_structure_track(fig, mv: MsaViz, profile: str, wrap_length: int):
    """Draw helix cylinders and strand arrows above each pyMSAviz MSA axis."""
    msa_axes = [ax for ax in fig.axes if tuple(round(v, 6) for v in ax.get_ylim()) == (0.0, float(mv.msa_count))]
    panel_count = (len(profile) + wrap_length - 1) // wrap_length
    if len(msa_axes) != panel_count:
        raise RuntimeError(
            f"Could not identify all wrapped MSA axes ({len(msa_axes)} found, "
            f"{panel_count} expected)."
        )

    for panel_index, ax in enumerate(msa_axes):
        start = panel_index * wrap_length
        end = min(start + wrap_length, len(profile))
        track_y = mv.msa_count + 0.18
        track_height = 0.55
        ax.set_ylim(0, mv.msa_count + 1.05)
        ax.text(start - 1, track_y + track_height / 2, "3BJX SS", ha="right", va="center", size=10)

        for left, right in runs(profile, HELIX_CODES, start, end):
            ax.add_patch(
                FancyBboxPatch(
                    (left, track_y), right - left, track_height,
                    boxstyle="round,pad=0.02,rounding_size=0.25",
                    facecolor="#d95f5f", edgecolor="#9e3030", linewidth=0.7,
                    clip_on=False, zorder=10,
                )
            )
        for left, right in runs(profile, STRAND_CODES, start, end):
            length = right - left
            head_length = min(1.4, length * 0.45)
            ax.add_patch(
                FancyArrow(
                    left, track_y + track_height / 2, length, 0,
                    width=track_height * 0.55, head_width=track_height,
                    head_length=head_length, length_includes_head=True,
                    facecolor="#4c78a8", edgecolor="#2f5680", linewidth=0.7,
                    clip_on=False, zorder=10,
                )
            )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("msa", nargs="?", type=Path, default=Path("dehi_first_20_aligned.fasta"))
    parser.add_argument("--pdb", type=Path, default=Path("3BJX.pdb"))
    parser.add_argument("--chain", default="A")
    parser.add_argument("--reference", default="3BJX")
    parser.add_argument(
        "--source", choices=("pdb", "dssp"), default="pdb",
        help="Use deposited PDB annotations (default) or recalculate with mkdssp",
    )
    parser.add_argument("--wrap-length", type=int, default=80)
    parser.add_argument("--output", type=Path, default=Path("dehi_3BJX_secondary_structure.png"))
    args = parser.parse_args()

    pdb_file = download_3bjx(args.pdb)
    alignment = AlignIO.read(args.msa, "fasta")
    reference = find_reference(alignment, args.reference)
    if args.source == "dssp":
        structure_sequence, structure_profile = dssp_sequence_and_profile(
            pdb_file, args.chain
        )
    else:
        structure_sequence, structure_profile = deposited_sequence_and_profile(
            pdb_file, args.chain
        )
    mapped_profile, identity = map_dssp_to_alignment(
        str(reference.seq), structure_sequence, structure_profile
    )
    if identity < 0.90:
        raise ValueError(
            f"3BJX MSA/DSSP sequence identity is only {identity:.1%}; "
            "check the reference row and chain."
        )

    mv = MsaViz(
        args.msa,
        wrap_length=args.wrap_length,
        show_count=True,
        show_grid=True,
        color_scheme="Clustal",
    )
    fig = mv.plotfig(dpi=150)
    add_secondary_structure_track(fig, mv, mapped_profile, args.wrap_length)
    fig.savefig(args.output, dpi=300, bbox_inches="tight", pad_inches=0.5)
    print(f"Mapped {args.source.upper()} annotations at {identity:.1%} sequence identity")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
