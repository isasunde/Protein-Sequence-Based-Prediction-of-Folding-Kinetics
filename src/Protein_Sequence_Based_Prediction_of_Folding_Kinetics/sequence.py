
import re
from functools import lru_cache

import pandas as pd
import requests


UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb"

MUTATION_PATTERN = re.compile(r"([A-Z])([1-9][0-9]*)([A-Z])")

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def parse_mutations(mutation_string):
    """Parse UniProt-style substitutions such as A34V or A34V,G56D."""

    if mutation_string is None or pd.isna(mutation_string):
        return []

    value = str(mutation_string).strip().strip("[]").strip()

    if not value or value.upper() == "WT":
        return []

    mutations = []

    for item in value.split(","):
        item = item.strip()
        match = MUTATION_PATTERN.fullmatch(item)

        if match is None:
            raise ValueError(f"Invalid mutation notation: {item!r}")

        wt, position, mutant = match.groups()

        if wt not in VALID_AMINO_ACIDS or mutant not in VALID_AMINO_ACIDS:
            raise ValueError(f"Invalid amino acid in mutation: {item!r}")

        mutations.append((wt, int(position), mutant))

    positions = [position for _, position, _ in mutations]

    if len(positions) != len(set(positions)):
        raise ValueError("Multiple substitutions at the same position.")

    return mutations


@lru_cache(maxsize=512)
def fetch_uniprot_sequence(accession):
    """Retrieve and cache a protein sequence from UniProt."""

    accession = str(accession).strip()

    if not accession or accession.lower() in {"nan", "none"}:
        raise ValueError("Missing UniProt accession.")

    url = f"{UNIPROT_FASTA_URL}/{accession}.fasta"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    lines = response.text.splitlines()

    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Unexpected FASTA response for {accession}")

    sequence = "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    ).upper()

    if not sequence:
        raise ValueError(f"Empty sequence for {accession}")

    invalid = set(sequence) - VALID_AMINO_ACIDS

    if invalid:
        raise ValueError(
            f"Unexpected sequence symbols for {accession}: {invalid}"
        )

    return sequence


def apply_mutations(sequence, mutation_string):
    """Apply substitutions using 1-based UniProt sequence positions."""

    mutations = parse_mutations(mutation_string)
    result = list(sequence)

    # Validate all substitutions against the same input sequence first.
    for wt, position, mutant in mutations:
        if position > len(sequence):
            raise ValueError(
                f"Position {position} exceeds sequence length {len(sequence)}"
            )

        observed = sequence[position - 1]

        if observed != wt:
            raise ValueError(
                f"Position {position}: expected {wt}, "
                f"but reference contains {observed}"
            )

    # Apply substitutions after validation.
    for wt, position, mutant in mutations:
        result[position - 1] = mutant

    return "".join(result)


def build_variant_sequences(
    reference_sequence,
    mutation_uniprot,
    background_mutation_uniprot=None,
):
    """
    Construct the experimental background and target variant sequences.

    Background substitutions are applied first. Target substitutions are
    then validated against that background.
    """

    background_sequence = apply_mutations(
        reference_sequence,
        background_mutation_uniprot,
    )

    mutant_sequence = apply_mutations(
        background_sequence,
        mutation_uniprot,
    )

    return {
        "sequence_reference": reference_sequence,
        "sequence_wildtype": background_sequence,
        "sequence_mutated": mutant_sequence,
    }


def extract_mutation_contexts(sequence, mutation_string, flank=15):
    """Return local sequence windows centred on all mutated positions."""

    starts = []
    ends = []

    for wt, position, mutant in parse_mutations(mutation_string):
        starts.append(max(1, position - flank))
        ends.append(min(len(sequence), position + flank))

    start = min(starts)
    end = max(ends)

    context = sequence[start - 1:end]

    return context
