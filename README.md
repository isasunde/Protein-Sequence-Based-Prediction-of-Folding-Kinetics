# Protein Sequence-Based Prediction of Folding Kinetics

## Overview

Protein folding kinetics describe how quickly a protein transitions from an unfolded state to its native, functional structure. Understanding the factors that influence folding rates can provide insight into protein stability, folding mechanisms, and the effects of mutations.

This project explores whether protein sequence information can be used to predict folding kinetics using experimental data from the K-Prot database. The workflow combines SQL-based data management, Python data processing, exploratory data analysis, and machine learning to investigate relationships between protein sequence, mutations, experimental conditions, and folding rates.

The long-term goal is to develop predictive models that can generalize to proteins not represented in the training data.

## Objectives

- **Data integration:** Retrieve and organize protein folding kinetics data from the K-Prot API and a relational MySQL database.
- **Data quality:** Validate records, handle missing values, standardize experimental metadata, and prepare a consistent modeling dataset.
- **Exploratory data analysis:** Investigate the distributions of folding rates and the influence of experimental conditions, measurement methods, and publications.
- **Sequence-based prediction:** Develop machine-learning models to predict folding kinetics from protein sequence information and mutation contexts.
- **Generalization:** Evaluate model performance on unseen proteins to assess whether the learned relationships extend beyond the training dataset.

## Data Source

The project uses experimental protein folding kinetics data from K-Prot, including protein identifiers, mutations, sequence information, kinetic measurements, and experimental metadata.

**Primary data source:** [K-Prot API](https://folding.biofold.org/k-pro/api)

The dataset includes measurements such as the logarithm of the folding rate extrapolated to water (`ln_kf_h2o`), alongside other kinetic and experimental variables where available.

## Methodology

### 1. Data acquisition and ETL

Data are retrieved from the K-Prot API and processed using a Python-based extract, transform, load (ETL) workflow.

The pipeline is designed to:
- Retrieve experimental records and relevant protein sequence information.
- Validate identifiers, mutations, and sequence references.
- Extract local sequence contexts around mutation sites based on UniProt IDs.
- Standardize data types and experimental metadata.
- Store structured records in a relational MySQL database.

### 2. Relational data management

A MySQL database organizes the data into related entities, including proteins, mutations, structures, publications, experiments, and kinetic measurements.

This structure separates protein-level information from experimental observations and supports reproducible data retrieval for downstream analysis.

### 3. Data quality and exploratory analysis

The analysis workflow retrieves a modeling dataset from SQL and assesses its completeness, consistency, and suitability for prediction.

Exploratory analysis focuses on:
- Folding-rate distributions and potential outliers.
- Missing data and coverage across proteins and mutations.
- Relationships between folding kinetics and experimental conditions, including temperature, pH, and measurement method.
- Differences between publications and potential confounding factors.
- The representation of wild-type and mutant sequences.

### 4. Feature engineering

Protein sequence information is transformed into features suitable for machine learning:
- Local sequence contexts surrounding mutation sites.
- Sequence-derived features describing amino acid substitutions.
- Protein language model representations, ESM embeddings.

These representations will be compared to investigate how different levels of sequence information contribute to predictive performance.

### 5. Model development and evaluation

Initial models will establish baseline performance before more complex sequence representations are introduced.

Evaluation focus on:
- Prediction error and agreement between observed and predicted folding rates.
- Comparisons between sequence-derived features and protein language model embeddings.
- Generalization to proteins excluded from model training.
- Potential data leakage arising from related mutations, repeated measurements, or shared publications.

To avoid artefacts from different data sources, different grouping approaches during train/test splitting is tested ensure that evaluation reflects the intended unseen-protein prediction task.

## Technology Stack

- **Programming:** Python
- **Data processing:** pandas, NumPy
- **Database:** MySQL, SQL
- **Data acquisition:** REST API, JSON
- **Visualization:** Matplotlib, Seaborn
- **Machine learning:** scikit-learn
- **Protein representations:** ESM protein language models
- **Development:** Git, Jupyter Notebook


## Getting Started

### Requirements

- Python 3.10 or later
- MySQL
- Access to the K-Prot dataset or a prepared local dataset

### Installation

Clone the repository and install the required Python packages:

```bash
git clone https://github.com/isasunde/Protein-Sequence-Based-Prediction-of-Folding-Kinetics.git
cd Protein-Sequence-Based-Prediction-of-Folding-Kinetics

pip install -r requirements.txt
```

Configure the database connection using environment variables. Create a local `.env` file based on `.env.example`:

```text
MYSQL_HOST=localhost
MYSQL_PORT=3306
MYSQL_USER=your_username
MYSQL_PASSWORD=your_password
MYSQL_DATABASE=kpro
```

### Running the analysis

Once the database is configured and populated, run the notebooks in numerical order, beginning with data quality assessment and exploratory analysis before proceeding to model development.

## Current Status

The project is being developed as a reproducible workflow from experimental data acquisition to protein-sequence-based prediction.

Current development priorities include:
- Validating the normalized database and measurement-level data integrity.
- Establishing a reliable SQL-to-Python analysis pipeline.
- Characterizing variation in folding kinetics across experimental conditions and publications.
- Building baseline prediction models before evaluating protein language model embeddings.

## Future Work

- Expand sequence feature engineering and compare alternative representations.
- Investigate the effects of mutations on folding-rate predictions.
- Evaluate generalization across proteins and publications.
- Quantify model uncertainty and identify limitations arising from heterogeneous experimental data.
- Improve reproducibility through automated data validation, testing, and documented analysis workflows.