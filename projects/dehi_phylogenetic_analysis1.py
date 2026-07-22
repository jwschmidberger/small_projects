"""Download protein sequences for this DEHI phylogenetic analysis.

This script can:

1. Download the amino-acid sequence for chain A of PDB entry 3BJX.
2. Download all UniProtKB amino-acid sequences matching Pfam family PF02627.
3. Download AhpD-family sequences using the specific InterPro family IPR004674.
4. Download DehI-like sequences using the specific InterPro family IPR019714.

The AhpD and DehI downloads deliberately use family assignments instead of
protein-name text, since names are not applied consistently across UniProtKB
records.
"""

from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen


PDB_ID = "3BJX"
CHAIN_ID = "A"
PFAM_ID = "PF02627"
AHPD_INTERPRO_ID = "IPR004674"
DEHI_INTERPRO_ID = "IPR019714"
PDB_OUTPUT_FASTA = Path(__file__).with_name(f"{PDB_ID}_{CHAIN_ID}.fasta")
UNIPROT_OUTPUT_FASTA = Path(__file__).with_name(f"uniprot_{PFAM_ID}.fasta")
AHPD_OUTPUT_FASTA = Path(__file__).with_name(
    f"uniprot_ahpd_{AHPD_INTERPRO_ID}.fasta"
)
DEHI_OUTPUT_FASTA = Path(__file__).with_name(
    f"uniprot_dehi_{DEHI_INTERPRO_ID}.fasta"
)
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


def download_uniprot_pfam_sequences(pfam_id: str = PFAM_ID) -> str:
    """Download all canonical UniProtKB sequences cross-referenced to a Pfam ID."""
    query = f"xref:pfam-{pfam_id.upper()}"
    url = f"{UNIPROT_API}?{urlencode({'query': query, 'format': 'fasta'})}"
    fasta = fetch_text(url)

    if not fasta.startswith(">"):
        raise RuntimeError(f"UniProt did not return FASTA data for {query}")

    UNIPROT_OUTPUT_FASTA.write_text(fasta)

    sequence_count = fasta.count("\n>")
    if fasta.startswith(">"):
        sequence_count += 1

    print(f"Wrote {sequence_count} UniProt sequences to {UNIPROT_OUTPUT_FASTA}")
    return fasta


def download_uniprot_ahpd_sequences(
    interpro_id: str = AHPD_INTERPRO_ID,
) -> str:
    """Download canonical UniProtKB sequences assigned to the AhpD family.

    IPR004674 is the AhpD family-level InterPro entry.  It is more specific than
    the CMD domain (PF02627/IPR003779), which also occurs outside AhpD proteins.
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
    2-haloacid dehalogenase DehI (also represented by Pfam PF10778).
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


def main() -> None:
    download_pdb_chain()
    download_uniprot_pfam_sequences()
    download_uniprot_ahpd_sequences()
    download_uniprot_dehi_sequences()



if __name__ == "__main__":
    main()
