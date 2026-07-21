"""Download the amino-acid sequence for chain A of PDB entry 3BJX.

The script uses the RCSB PDB REST API and writes the sequence as FASTA.
"""

from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


PDB_ID = "3BJX"
CHAIN_ID = "A"
OUTPUT_FASTA = Path(__file__).with_name(f"{PDB_ID}_{CHAIN_ID}.fasta")
RCSB_API = "https://data.rcsb.org/rest/v1/core"


def fetch_json(url: str) -> dict:
    """Fetch JSON from an RCSB API endpoint."""
    try:
        with urlopen(url, timeout=30) as response:
            return json.load(response)
    except HTTPError as exc:
        raise RuntimeError(f"RCSB returned HTTP {exc.code} for {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach RCSB PDB at {url}: {exc.reason}") from exc


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


def main() -> None:
    sequence = sequence_for_chain(PDB_ID, CHAIN_ID)
    write_fasta(sequence, OUTPUT_FASTA, PDB_ID, CHAIN_ID)
    print(f"Wrote {len(sequence)} amino acids to {OUTPUT_FASTA}")


if __name__ == "__main__":
    main()
