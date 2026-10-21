import re
import requests
import pandas as pd
from datetime import datetime

from .sequence import (
    fetch_uniprot_sequence,
    build_variant_sequences,
    extract_mutation_contexts,
)

# Prepare regex for pattern recognition of set columns
UNIPROT_REGEX = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)

MUTATION_REGEX = re.compile(
    r'(?:[A-Z][1-9][0-9]*[A-Z])(?:,[A-Z][1-9][0-9]*[A-Z])*|WT'
)

single_mutation_regex = r'[A-Z][1-9][0-9]*[A-Z]'
double_mutation_regex = rf'{single_mutation_regex},{single_mutation_regex}'

NUMERIC_COLUMNS = ["temperature", "ph", "ln_kf", "ln_ku", "delta_ln_kf", "delta_ln_ku", "rsa"]

# Expected columns of K-pro website
EXP_COLUMNS = ['PROTEIN', 'SOURCE', 'LENGTH', 'UniProt', 'PFAM',
               'CATH', 'EC_NUMBER', 'PDB_wild', 'MUTATED_CHAIN',
               'MUTATION_UNIPROT', 'MUTATION_PDB', 'SEC_STR',
               'RSA', 'T', 'pH', 'BUFFER_NAME', 'BUFFER_CONC',
               'ION_NAME', 'ADDITIVES', 'MEASURE_KIN', 'METHOD_KIN',
               'ln(kf)_H2O', 'dln(kf)_H2O', 'ln(ku)_H2O', 'dln(ku)_H2O',
               'ln(kf)_DEN', 'dln(kf)_DEN', 'ln(ku)_DEN', 'dln(ku)_DEN',
               'dGKIN_H2O', 'ddGKIN_H2O', 'mf', 'mu', 'CM_KIN', 'beta-T',
               'PhiF_H2O', 'STATE', 'REVERSIBILITY', 'AUTHOR_KIN',
               'REFERENCE_KIN', 'YEAR_KIN', 'PMID_KIN', 'REVIEW_DATE',
               'BACKGROUND_MUT_UNIPROT', 'BACKGROUND_MUT_PDB'
               ]

def prepare_columns(data):
    """Add any missing expected columns with None values."""
    # Find columns not present
    missing_columns = [
        column for column in EXP_COLUMNS
        if column not in data.columns
    ]
    # Add missing column to dataframe
    if missing_columns:
        data = data.copy()
        data[missing_columns] = None

    return data

def fill_protein_identifiers(data):
    """Add protein name from uniprot ID or reverse, if one is missing"""
    data = data.copy()

    # Normalize empty strings to missing values
    for column in ["PROTEIN", "UniProt"]:
        data[column] = data[column].replace(r"^\s*$", pd.NA, regex=True)

    # Fill missing UniProt IDs using the protein name
    uniprot_by_protein = (
        data.dropna(subset=["PROTEIN", "UniProt"])
        .drop_duplicates("PROTEIN")
        .set_index("PROTEIN")["UniProt"]
    )

    missing_uniprot = data["UniProt"].isna()
    data.loc[missing_uniprot, "UniProt"] = (
        data.loc[missing_uniprot, "PROTEIN"].map(uniprot_by_protein)
    )

    # Fill missing protein names using the UniProt ID
    protein_by_uniprot = (
        data.dropna(subset=["UniProt", "PROTEIN"])
        .drop_duplicates("UniProt")
        .set_index("UniProt")["PROTEIN"]
    )

    missing_protein = data["PROTEIN"].isna()
    data.loc[missing_protein, "PROTEIN"] = (
        data.loc[missing_protein, "UniProt"].map(protein_by_uniprot)
    )

    # Use categorical dtype if desired
    data["UniProt"] = data["UniProt"].astype("category")

    return data

def align_dtype(data):
    """Set d-type to the expected for all columns"""

    categorical_columns = ['PROTEIN', 'SOURCE', 'EC_NUMBER', 
                           'PDB_wild', 'MUTATED_CHAIN', 'SEC_STR',
                           'STATE', 'REVERSIBILITY', 'AUTHOR_KIN',
                           'REFERENCE_KIN',]

    integer_columns = ['LENGTH', 'YEAR_KIN', 'PMID_KIN']

    float_columns = ['T', 'pH', 'BUFFER_CONC', 
                     'ln(kf)_H2O', 'ln(ku)_H2O', 'dln(kf)_H2O', 
                     'dln(ku)_H2O', 'ln(kf)_DEN', 'ln(ku)_DEN', 
                     'dln(kf)_DEN', 'dln(ku)_DEN', 'dGKIN_H2O', 
                     'ddGKIN_H2O', 'mf', 'mu', 'CM_KIN', 'beta-T'
                     ]

    string_columns = ['PFAM', 'CATH', 'MUTATION_UNIPROT', 
                      'MUTATION_PDB', 'BUFFER_NAME', 'ION_NAME', 
                      'ADDITIVES', 'MEASURE_KIN', 'METHOD_KIN',
                      'BACKGROUND_MUT_UNIPROT', 'BACKGROUND_MUT_PDB',
                      'RSA']

    # Input UNIPROT ID or Protein name from other entries of the same protein when missing
    data = fill_protein_identifiers(data)
    data['REVERSIBILITY'] = data['REVERSIBILITY'].str.title()

    # Clean known issue value (column contain >1 values, which are set to NaN)
    data['PhiF_H2O'] = pd.to_numeric(data['PhiF_H2O'], errors = 'coerce')

    try:
        for column in categorical_columns:
            data[column] = data[column].fillna('NA')
            data[column] = data[column].replace(["", "NA", "N/A", "None", 'NaN'], 'NA')
            data[column] = data[column].astype('category')

        for column in integer_columns:
            data[column] = data[column].astype(int)
            data[column] = data[column].replace(["", "NA", "N/A", "None"], pd.NA)

        for column in float_columns:
            data[column] = data[column].astype(float)
            data[column] = data[column].replace(["", "NA", "N/A", "None"], pd.NA)
    
        for column in string_columns:
            data[column] = data[column].fillna('NA')
            data[column] = data[column].astype(str)
            data[column] = data[column].replace(["", "NA", "N/A", "None"], 'NA')
        
    except ValueError:
        print(f'ValueError upon conversion of {column}')

    data['REVIEW_DATE'] = pd.to_datetime(data['REVIEW_DATE'])

    # Remove parenthesis
    data['PFAM'] = data['PFAM'].str.strip('[').str.strip(']')
    data['CATH'] = data['CATH'].str.strip('[').str.strip(']')
    data['BACKGROUND_MUT_UNIPROT'] = data['BACKGROUND_MUT_UNIPROT'].str.strip('[').str.strip(']')
    data['BACKGROUND_MUT_PDB'] = data['BACKGROUND_MUT_PDB'].str.strip('[').str.strip(']')

    # Replace inconsistent buffer names
    BUFFER_MAPPING = {
        "Na-Acetate": "Sodium Acetate",
        "NaAcetate": "Sodium Acetate",
        "Na-acetate": "Sodium Acetate",
        "Na Acetate": "Sodium Acetate",
    }
    data["BUFFER_NAME"] = data["BUFFER_NAME"].replace(BUFFER_MAPPING)

    MEASURE_KIN_MAPPING = {
        "stopped-flow": "Stopped-flow",
        "T-jump": "T-jump",
        "T-jump or continuous-flow": "T-jump or continuous-flow",
        "T-jump or stopped-flow": "T-jump or stopped-flow",
        "stopped-flow, or T-Jump Fluorescence Sectroscopy":
            "Stopped-flow or T-jump fluorescence spectroscopy",
        "stopped-flow or T-Jump Fluorescence Sectroscopy":
            "Stopped-flow or T-jump fluorescence spectroscopy",
        "manual mixing": "Manual mixing",
        "continuous-flow": "Continuous-flow",
    }
    data["MEASURE_KIN"] = data["MEASURE_KIN"].str.strip().replace(MEASURE_KIN_MAPPING)

    
    METHOD_KIN_MAPPING = {
        "GuHCl": "Guanidine hydrochloride",
        "GdmCl": "Guanidine hydrochloride",
        "GndHCl": "Guanidine hydrochloride",
        "GdnHCl": "Guanidine hydrochloride",

        "urea": "Urea",
        "Urea": "Urea",

        "thermal": "Thermal denaturation",
        "thermal denaturation": "Thermal denaturation",
        "heat": "Thermal denaturation",

        "thermal or GuHCl":
            "Thermal denaturation or guanidine hydrochloride",
        "alkaline unfolding (pH 12.4)":
            "Alkaline unfolding (pH 12.4)",
    }
    data["METHOD_KIN"] = data["METHOD_KIN"].str.strip().replace(METHOD_KIN_MAPPING)
    
    return data

def add_sequences(
    data: pd.DataFrame,
    flank: int = 100,
) -> pd.DataFrame:
    """
    Add local amino-acid sequence contexts surrounding mutations.

    Reference and mutant sequences are constructed temporarily for
    validation, but only the local mutation contexts are retained.

    Parameters
    ----------
    data : pd.DataFrame
        K-Pro records containing UniProt, MUTATION_UNIPROT and,
        optionally, BACKGROUND_MUT_UNIPROT.
    flank : int, default=100
        Number of amino acids to include on each side of each mutation.

    Returns
    -------
    pd.DataFrame
        Copy of the input DataFrame with mutation_contexts and
        sequence_status columns added.
    """

    data = data.copy()
    data["reference_context"] = None
    data["mutation_context"] = None
    data["sequence_status"] = None

    # Cache reference sequences so each accession is retrieved only once.
    sequence_cache = {}

    for accession in data["UniProt"].dropna().unique():
        try:
            sequence_cache[accession] = fetch_uniprot_sequence(accession)
        except (requests.RequestException, ValueError) as exc:
            sequence_cache[accession] = exc

    for idx, row in data.iterrows():
        accession = row["UniProt"]

        if pd.isna(accession):
            data.at[idx, "sequence_status"] = "missing_accession"
            continue

        reference = sequence_cache.get(accession)

        if isinstance(reference, Exception):
            data.at[idx, "sequence_status"] = (
                f"reference_retrieval_failed: {reference}"
            )
            continue

        if reference is None:
            data.at[idx, "sequence_status"] = "reference_not_found"
            continue

        try:
            sequences = build_variant_sequences(
                reference_sequence=reference,
                mutation_uniprot=row["MUTATION_UNIPROT"],
                background_mutation_uniprot=row.get(
                    "BACKGROUND_MUT_UNIPROT"
                ),
            )

            reference_context = extract_mutation_contexts(
                sequence=sequences["sequence_wildtype"],
                mutation_string=row["MUTATION_UNIPROT"],
                flank=flank,
            )
            if row['MUTATION_UNIPROT'] == 'WT':
                mutation_context = None
                status = 'No mutation'
                
            else:
                mutation_context = extract_mutation_contexts(
                    sequence=sequences["sequence_mutated"],
                    mutation_string=row["MUTATION_UNIPROT"],
                    flank=flank,
                )
                status = "validated"
                
            data.at[idx, "reference_context"] = reference_context
            data.at[idx, "mutation_context"] = mutation_context
            data.at[idx, "sequence_status"] = status

        except ValueError as exc:
            data.at[idx, "sequence_status"] = (
                f"mutation_validation_failed: {exc}"
            )

    return data


def clean(data: pd.DataFrame) -> pd.DataFrame:
    """Clean and standardize K-Pro data."""

    data = data.copy()

    # Align data types
    data = align_dtype(data)

    # Remove complete duplicates
    data.drop_duplicates(subset = ['UniProt', 'SOURCE', 'BUFFER_NAME', 'BUFFER_CONC', 'MUTATION_PDB', 'PDB_wild', 'T', 'pH', 'METHOD_KIN', 'MEASURE_KIN', 'STATE', 'PMID_KIN'],inplace = True)

    # Strip whitespace from string columns
    string_columns = data.select_dtypes(include=['object', 'str']).columns
    data[string_columns] = data[string_columns].apply(
        lambda column: column.str.strip()
    )

    # Fill missing mutation identifiers
    data['MUTATION_PDB'] = data['MUTATION_PDB'].fillna(data['MUTATION_UNIPROT'])
    data['MUTATION_UNIPROT'] = data['MUTATION_UNIPROT'].fillna(data['MUTATION_PDB'])
    return data

def transform(data: pd.DataFrame) -> pd.DataFrame:
    """Transform cleaned K-Pro records into a standardized representation."""

    data = data.copy()

    data = data.rename(columns={
        "PROTEIN": "protein_name",
        "SOURCE": "protein_source",
        "LENGTH": "protein_length",
        "UniProt": "uniprot_id",
        "PDB_wild": "pdb_id",
        "MUTATED_CHAIN": "mutated_chain",
        "PFAM": "pfam",
        "CATH": "cath",
        "EC_NUMBER": "ec_number",
        "MUTATION_UNIPROT": "mutation_uniprot",
        "MUTATION_PDB": "mutation_pdb",
        "T": "temperature",
        "pH": "ph",
        "BUFFER_NAME": "buffer_name",
        "BUFFER_CONC": "buffer_conc",
        "ION_NAME": 'ion_name',
        "MEASURE_KIN": "kinetic_measurement",
        "METHOD_KIN": "kinetic_method",
        "ln(kf)_H2O": "ln_kf_h2o",
        "ln(ku)_H2O": "ln_ku_h2o",
        "dGKIN_H2O": "dg_kin_h2o",
        "ln(kf)_DEN": 'ln_kf_den',
        "ln(ku)_DEN": 'ln_ku_den',
        "mf": "mf",
        "mu": "mu",
        "AUTHOR_KIN": "authors",
        "REFERENCE_KIN": "ref",
        "YEAR_KIN": "publication_year",
        "PMID_KIN": "pmid",
        "STATE": "state",
        "REVERSIBILITY": "reversibility",
        "REVIEW_DATE": "review_date",
        "ADDITIVES": "additives",
        "beta-T": "beta_t",
        "SEC_STR": "secondary_structure",
        "RSA": "rsa",
        "CM_KIN": "cm",
        "dln(kf)_H2O": "dln_kf_h2o",
        "dln(ku)_H2O": "dln_ku_h2o", 
        "ddGKIN_H2O": "ddgkin_h2o",
        "dln(kf)_DEN": 'dln_kf_den',
        "dln(ku)_DEN": 'dln_ku_den',
        "PhiF_H2O": "phif_h2o", 
        'BACKGROUND_MUT_UNIPROT':'background_mutation_uniprot',
        'BACKGROUND_MUT_PDB':'background_mutation_pdb'
    })

    for column in NUMERIC_COLUMNS:
            if column in data.columns:
                data[column] = pd.to_numeric(
                    data[column],
                    errors="coerce"
                )

    return data

def validate(data: pd.DataFrame) -> dict:
    """Check essential requirements before loading records into SQL."""

    if data.empty:
        raise ValueError("ETL produced no records.")

    required = ["uniprot_id", "reference_context"]
    missing = set(required) - set(data.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    missing_target = data["ln_kf_h2o"].isna()
    missing_id = data["uniprot_id"].isna()

    sequence_ok = data["sequence_status"].eq("validated")
    summary = {
        "records": len(data),
        "missing_targets": int(missing_target.sum()),
        "missing_identifiers": int(missing_id.sum()),
        "valid_sequences": int(sequence_ok.sum()),
        "invalid_sequences": int((~sequence_ok).sum()),
    }

    if missing_id.any():
        raise ValueError(
            f"{missing_id.sum()} records have no UniProt identifier."
        )

    return summary

def quality_assurance(data: pd.DataFrame):
    print('='*30)
    print('K-pro Data Quality Report')
    print('='*30)

    print('Records: '+ str(data['UniProt'].count()))
    print('Mutations: '+ str(data['MUTATION_PDB'].nunique()))
    print('Publications: '+ str(data['PMID_KIN'].nunique()))

    print('\nMissing values')
    print('-'*10)

    na_columns = data.isna().sum()
    for i, no in na_columns.items():
        if no != 0:
            print(i + ': ' + str(no))

    print('\nMutation types')
    print('-'*10)
    print('WT: ' + str((data['MUTATION_PDB'] == 'WT').sum()))
    print('Single: ' + str(data['MUTATION_PDB'].str.fullmatch(single_mutation_regex).sum()))
    print('Double: ' + str(data['MUTATION_PDB'].str.fullmatch(double_mutation_regex).sum()))

    print('\nExperimental Conditions')
    print('-'*10)
    print('Temperature: '+ str(data['T'].nunique()) + ' unique values')
    print('pH: '+ str(data['pH'].nunique()) + ' unique values')
    print('Methods: '+ str(data['METHOD_KIN'].nunique()) + ' unique values')
    print('Buffer: '+ str(data['BUFFER_NAME'].nunique()) + ' unique values')
    print('States: '+ str(data['STATE'].nunique()) + ' unique values')

    print('\nDuplicates')
    print('Exact duplicates: ' + str(data.drop(columns = ['CATH', 'PFAM']).duplicated().sum()))
    print('Potential duplicates: ' + str(data.duplicated(subset = ['UniProt', 'SOURCE', 'BUFFER_NAME', 'BUFFER_CONC', 'MUTATION_PDB', 'PDB_wild', 'T', 'pH', 'METHOD_KIN', 'MEASURE_KIN', 'STATE', 'PMID_KIN']).sum()))
