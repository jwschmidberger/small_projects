"""Collect sequences for a CMD/AhpD/DehI evolutionary analysis.

Research question
-----------------
Test the hypothesis that DehI-like proteins arose through duplication and
fusion of a smaller ancestral unit related to proteins in the CMD and/or AhpD
groups.  The three sequence collections produced here are the input data for
subsequent domain-boundary analysis, multiple-sequence alignment, and
phylogenetic reconstruction; this downloader does not itself test the
hypothesis.

Sequence collections
--------------------

1. CMD-domain proteins: UniProtKB records matching Pfam PF02627.
2. AhpD-family proteins: UniProtKB records matching InterPro IPR004674.
3. DehI-like proteins: UniProtKB records matching InterPro IPR019714.
4. The sequence of chain A from the DehI reference structure PDB 3BJX.

Each of the three family collections is clustered independently at a chosen
amino-acid identity (90% by default) with MMseqs2.  Independent clustering
preserves the family label needed for later phylogenetic comparisons.  At 90%,
the representative sequences are written to ``*_nr90.fasta`` files, while
``*_nr90_clusters.tsv`` files record the representative/member mapping.

Database family assignments are used instead of protein-name searches because
annotation names vary.  Full-length sequences are retained so a downstream
analysis can compare the N- and C-terminal regions of DehI separately with the
smaller CMD and AhpD proteins, as required to evaluate a duplication/fusion
model.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen


PDB_ID = "3BJX"
CHAIN_ID = "A"
# Family identifiers defining the three comparison groups.  PF02627 is the CMD
# domain; IPR004674 and IPR019714 are the specific AhpD and DehI families.
CMD_PFAM_ID = "PF02627"
AHPD_INTERPRO_ID = "IPR004674"
DEHI_INTERPRO_ID = "IPR019714"
# Retained as a compatibility name for callers of the original CMD downloader.
PFAM_ID = CMD_PFAM_ID

# 3BJX is an experimentally determined DehI structure and provides a useful
# reference sequence for locating a possible internal duplication boundary.
RESULTS_DIR = Path(__file__).with_name("results")
PDB_OUTPUT_FASTA = RESULTS_DIR / f"{PDB_ID}_{CHAIN_ID}.fasta"
CMD_OUTPUT_FASTA = RESULTS_DIR / f"uniprot_{CMD_PFAM_ID}.fasta"
# Retained as a compatibility name for code using the original output constant.
UNIPROT_OUTPUT_FASTA = CMD_OUTPUT_FASTA
AHPD_OUTPUT_FASTA = RESULTS_DIR / f"uniprot_ahpd_{AHPD_INTERPRO_ID}.fasta"
DEHI_OUTPUT_FASTA = RESULTS_DIR / f"uniprot_dehi_{DEHI_INTERPRO_ID}.fasta"
FAMILY_FASTA_FILES = (CMD_OUTPUT_FASTA, AHPD_OUTPUT_FASTA, DEHI_OUTPUT_FASTA)
CLUSTER_IDENTITY = 0.90
# Bidirectional coverage over 80% of the longer sequence keeps fragments from
# collapsing full-length proteins and helps preserve multidomain architecture.
MIN_ALIGNMENT_COVERAGE = 0.80
RCSB_API = "https://data.rcsb.org/rest/v1/core"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/stream"


def fetch_json(url: str) -> dict:
    """Fetch JSON from an API endpoint."""
    try:
        with urlopen(url, timeout=30) as response:
            return json.load(response)
    except HTTPError as exc:
        raise RuntimeError(f"RCSB returned HTTP {exc.code} for {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach RCSB PDB at {url}: {exc.reason}") from exc


def fetch_text(url: str) -> str:
    """Fetch plain text from an API endpoint."""
    try:
        with urlopen(url, timeout=120) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        details = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"API returned HTTP {exc.code} for {url}\n{details}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach API at {url}: {exc.reason}") from exc


def sequence_for_chain(pdb_id: str, chain_id: str) -> str:
    """Return the canonical one-letter amino-acid sequence for a PDB chain."""
    pdb_id = pdb_id.upper()

    instance_url = f"{RCSB_API}/polymer_entity_instance/{pdb_id}/{chain_id}"
    instance = fetch_json(instance_url)

    identifiers = instance["rcsb_polymer_entity_instance_container_identifiers"]
    entity_id = identifiers["entity_id"]

    entity_url = f"{RCSB_API}/polymer_entity/{pdb_id}/{entity_id}"
    entity = fetch_json(entity_url)

    sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
    return "".join(sequence.split())


def write_fasta(sequence: str, output_path: Path, pdb_id: str, chain_id: str) -> None:
    """Write a wrapped FASTA file."""
    header = f">{pdb_id.upper()}_{chain_id} chain {chain_id} amino-acid sequence"
    wrapped = "\n".join(
        sequence[position : position + 80] for position in range(0, len(sequence), 80)
    )
    output_path.write_text(f"{header}\n{wrapped}\n")


def download_pdb_chain() -> str:
    sequence = sequence_for_chain(PDB_ID, CHAIN_ID)
    write_fasta(sequence, PDB_OUTPUT_FASTA, PDB_ID, CHAIN_ID)
    print(f"Wrote {len(sequence)} amino acids to {PDB_OUTPUT_FASTA}")
    return sequence


def download_uniprot_cmd_sequences(pfam_id: str = CMD_PFAM_ID) -> str:
    """Download the CMD comparison group using its Pfam assignment.

    CMD proteins represent one set of smaller structural relatives against
    which the two putative halves of DehI can be compared downstream.
    """
    query = f"xref:pfam-{pfam_id.upper()}"
    url = f"{UNIPROT_API}?{urlencode({'query': query, 'format': 'fasta'})}"
    fasta = fetch_text(url)

    if not fasta.startswith(">"):
        raise RuntimeError(f"UniProt did not return FASTA data for {query}")

    CMD_OUTPUT_FASTA.write_text(fasta)

    sequence_count = fasta.count("\n>")
    if fasta.startswith(">"):
        sequence_count += 1

    print(f"Wrote {sequence_count} CMD-domain sequences to {CMD_OUTPUT_FASTA}")
    return fasta


def download_uniprot_pfam_sequences(pfam_id: str = PFAM_ID) -> str:
    """Backward-compatible name for :func:`download_uniprot_cmd_sequences`."""
    return download_uniprot_cmd_sequences(pfam_id)


def download_uniprot_ahpd_sequences(
    interpro_id: str = AHPD_INTERPRO_ID,
) -> str:
    """Download canonical UniProtKB sequences assigned to the AhpD family.

    IPR004674 is the AhpD family-level InterPro entry.  It is more specific than
    the CMD domain (PF02627/IPR003779), which also occurs outside AhpD proteins.
    AhpD proteins provide the second group of smaller candidate relatives for
    comparison with the putative repeated regions of DehI.
    """
    interpro_id = interpro_id.upper()
    query = f"xref:interpro-{interpro_id}"
    url = f"{UNIPROT_API}?{urlencode({'query': query, 'format': 'fasta'})}"
    fasta = fetch_text(url)

    if not fasta.startswith(">"):
        raise RuntimeError(f"UniProt did not return FASTA data for {query}")

    AHPD_OUTPUT_FASTA.write_text(fasta)
    sequence_count = fasta.count("\n>") + 1
    print(f"Wrote {sequence_count} AhpD-family sequences to {AHPD_OUTPUT_FASTA}")
    return fasta


def download_uniprot_dehi_sequences(
    interpro_id: str = DEHI_INTERPRO_ID,
) -> str:
    """Download canonical UniProtKB sequences assigned to the DehI family.

    IPR019714 is the family-level InterPro entry for configuration-inverting
    2-haloacid dehalogenase DehI (also represented by Pfam PF10778).  These
    full-length proteins are the proposed duplication/fusion products; a later
    analysis should infer or define their internal boundary before comparing
    each region with CMD and AhpD sequences.
    """
    interpro_id = interpro_id.upper()
    query = f"xref:interpro-{interpro_id}"
    url = f"{UNIPROT_API}?{urlencode({'query': query, 'format': 'fasta'})}"
    fasta = fetch_text(url)

    if not fasta.startswith(">"):
        raise RuntimeError(f"UniProt did not return FASTA data for {query}")

    DEHI_OUTPUT_FASTA.write_text(fasta)
    sequence_count = fasta.count("\n>") + 1
    print(f"Wrote {sequence_count} DehI-family sequences to {DEHI_OUTPUT_FASTA}")
    return fasta


def count_fasta_sequences(path: Path) -> int:
    """Count FASTA records without loading the complete collection into memory."""
    with path.open() as fasta_file:
        return sum(line.startswith(">") for line in fasta_file)


def identity_label(identity: float) -> str:
    """Return an output label such as ``nr90`` for an identity fraction."""
    if not 0.5 <= identity <= 1.0:
        raise ValueError("Clustering identity must be between 0.5 and 1.0")
    return f"nr{identity * 100:g}"


def cluster_sequence_collection(
    input_path: Path, identity: float = CLUSTER_IDENTITY
) -> Path:
    """Cluster one protein collection at the requested identity threshold.

    MMseqs2 ``easy-linclust`` scales linearly with collection size.  Exact
    alignment identity is used instead of MMseqs2's default score-derived
    estimate, and 80% bidirectional coverage protects full-length architecture.
    """
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        raise RuntimeError(
            "MMseqs2 is required for 90% clustering but was not found on PATH. "
            "Install it with 'brew install mmseqs2' or your system package manager."
        )
    if not input_path.is_file():
        raise FileNotFoundError(f"Sequence collection not found: {input_path}")

    label = identity_label(identity)
    output_path = input_path.with_name(f"{input_path.stem}_{label}.fasta")
    membership_path = input_path.with_name(
        f"{input_path.stem}_{label}_clusters.tsv"
    )
    with tempfile.TemporaryDirectory(
        prefix=f"{input_path.stem}_mmseqs_", dir=input_path.parent
    ) as temporary_directory:
        temporary_path = Path(temporary_directory)
        result_prefix = temporary_path / "clusters"
        command = [
            mmseqs,
            "easy-linclust",
            str(input_path),
            str(result_prefix),
            str(temporary_path / "work"),
            "--min-seq-id",
            str(identity),
            "-c",
            str(MIN_ALIGNMENT_COVERAGE),
            "--cov-mode",
            "0",
            "--alignment-mode",
            "3",
            "-v",
            "1",
        ]
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"MMseqs2 failed while clustering {input_path}"
            ) from exc

        shutil.copyfile(f"{result_prefix}_rep_seq.fasta", output_path)
        shutil.copyfile(f"{result_prefix}_cluster.tsv", membership_path)

    input_count = count_fasta_sequences(input_path)
    representative_count = count_fasta_sequences(output_path)
    print(
        f"Reduced {input_path.name} from {input_count} to "
        f"{representative_count} representatives: {output_path}"
    )
    return output_path


def cluster_family_collections(
    identity: float = CLUSTER_IDENTITY,
) -> list[Path]:
    """Cluster the CMD, AhpD, and DehI groups independently."""
    return [
        cluster_sequence_collection(path, identity) for path in FAMILY_FASTA_FILES
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cluster-only",
        action="store_true",
        help="cluster the existing family FASTA files without downloading them again",
    )
    parser.add_argument(
        "--cluster-identity",
        type=float,
        default=CLUSTER_IDENTITY,
        metavar="FRACTION",
        help="minimum sequence identity as a fraction from 0.5 to 1.0 (default: 0.90)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        identity_label(args.cluster_identity)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if not args.cluster_only:
        download_pdb_chain()
        download_uniprot_cmd_sequences()
        download_uniprot_ahpd_sequences()
        download_uniprot_dehi_sequences()
    cluster_family_collections(args.cluster_identity)



if __name__ == "__main__":
    main()
