import re
import pandas as pd
from datetime import datetime

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

    print(f"Added columns:    {missing_columns}")

    return data

def fill_protein_identifiers(data):
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
                     'ddGKIN_H2O', 'mf', 'mu', 'CM_KIN', 'beta-T',
                     'PhiF_H2O',
                     ]

    string_columns = ['PFAM', 'CATH', 'MUTATION_UNIPROT', 
                      'MUTATION_PDB', 'BUFFER_NAME', 'ION_NAME', 
                      'ADDITIVES', 'MEASURE_KIN', 'METHOD_KIN',
                      'BACKGROUND_MUT_UNIPROT', 'BACKGROUND_MUT_PDB',
                      'RSA']

    # Input UNIPROT ID or Protein name from other entries of the same protein when missing
    data = fill_protein_identifiers(data)
    data['REVERSIBILITY'] = data['REVERSIBILITY'].str.title()

    try:
        for column in categorical_columns:
            data[column] = data[column].astype('category')

        for column in integer_columns:
            data[column] = data[column].astype(int)

        for column in float_columns:
            data[column] = data[column].astype(float)
    
        for column in string_columns:
            data[column] = data[column].astype(str)
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


def validate_record(entry):
    """Validate a single K-Pro record and return any validation errors."""

    errors = []

    # Protein
    # ---------------------------------------
    protein = entry.get("PROTEIN")

    if protein is None or not str(protein).strip():
        errors.append("Missing protein")

    # UniProt ID
    # ---------------------------------------
    uniprot = entry.get("UniProt")

    if not isinstance(uniprot, str) or not UNIPROT_REGEX.fullmatch(uniprot):
        errors.append("Invalid UniProt")

    # Temperature
    # ---------------------------------------
    temperature = entry.get("T")

    try:
        temperature = float(temperature)

        if not -273.15 <= temperature <= 100:
            errors.append("Temperature outside valid range")

    except (TypeError, ValueError):
        errors.append("Invalid temperature")

    # pH
    # ---------------------------------------
    ph = entry.get("pH")

    try:
        ph = float(ph)

        if not 0 <= ph <= 14:
            errors.append("pH outside valid range")

    except (TypeError, ValueError):
        errors.append("Invalid pH")

    # RSA
    # ---------------------------------------
    rsa = entry.get("RSA")
    
    try:
        rsa = float(rsa)

        if not 0 <= rsa <= 100:
            errors.append("RSA outside valid range")
    
    except (TypeError, ValueError):
        errors.append("Invalid RSA")


    # Mutations
    # ---------------------------------------
    mutation_pdb = entry.get("MUTATION_PDB")

    if not isinstance(mutation_pdb, str) or not MUTATION_REGEX.fullmatch(mutation_pdb):
        errors.append("Invalid PDB mutation")

    mutation_uniprot = entry.get("MUTATION_UNIPROT")

    if not isinstance(mutation_uniprot, str) or not MUTATION_REGEX.fullmatch(mutation_uniprot):
        errors.append("Invalid UniProt mutation")

    return errors

def validation_report(data: pd.DataFrame) -> pd.DataFrame:
    """Generate a validation report for the extracted K-Pro data."""

    validation_results = data.apply(
        validate_record, axis=1)

    rules = [
        "Missing protein",
        "Invalid UniProt",
        "Invalid temperature",
        "Temperature outside valid range",
        "Invalid pH",
        "pH outside valid range",
        "Invalid RSA",
        "RSA outside valid range",
        "Invalid PDB mutation",
        "Invalid UniProt mutation"
    ]

    report = []

    total_records = len(data)

    for rule in rules:
        failed = validation_results.apply(
            lambda errors: rule in errors
            ).sum()

        report.append({
            "Rule": rule,
            "Failed": failed,
            "Passed": total_records - failed,
            "Failure %": 100 * failed / total_records
        })

    return pd.DataFrame(report)

def clean(data: pd.DataFrame) -> pd.DataFrame:
    """Clean and standardize K-Pro data without changing its meaning."""

    data = data.copy()

    # Remove complete duplicates
    data.drop_duplicates(subset = ['UniProt', 'SOURCE', 'BUFFER_NAME', 'BUFFER_CONC', 'MUTATION_PDB', 'PDB_wild', 'T', 'pH', 'METHOD_KIN', 'MEASURE_KIN', 'STATE', 'PMID_KIN'],inplace = True)
    
    # Standardize missing values
    data = data.replace(["", "NA", "N/A", "None"], pd.NA)

    # Strip whitespace from string columns
    string_columns = data.select_dtypes(include=['object', 'str']).columns
    data[string_columns] = data[string_columns].apply(
        lambda column: column.str.strip()
    )

    # Convert numeric columns
    numeric_columns = ["T", "pH", "RSA", "LENGTH", "ln(kf)_H2O", "ln(ku)_H2O",
                       "dGKIN_H2O", "mf", "mu", "beta-T",'PhiF_H2O']
    
    data[numeric_columns] = data[numeric_columns].apply(
        pd.to_numeric,
        errors="coerce"
    )

    return data

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
