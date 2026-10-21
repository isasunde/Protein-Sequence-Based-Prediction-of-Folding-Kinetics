# Create database
DROP DATABASE kpro;
CREATE DATABASE IF NOT EXISTS kpro;

# Use created database
USE kpro;

CREATE TABLE IF NOT EXISTS placeholder(
    record INT PRIMARY KEY AUTO_INCREMENT,
    protein_name VARCHAR(200),
    protein_source VARCHAR(100),
    protein_length DECIMAL(5, 2),
    uniprot_id VARCHAR(10),
    pdb_id VARCHAR(12),
    mutated_chain VARCHAR(2),
    pfam VARCHAR(50),
    cath VARCHAR(100),
    ec_number VARCHAR(10),
    mutation_uniprot VARCHAR(50),
    mutation_pdb VARCHAR(50),
    temperature DECIMAL(5, 2),
    ph DECIMAL(5, 2),
    buffer_name VARCHAR(50),
    buffer_conc VARCHAR(50),
    ion_name VARCHAR(50),
    kinetic_measurement VARCHAR(50),
    kinetic_method VARCHAR(50),
    ln_kf_h2o DECIMAL(5, 2),
    ln_ku_h2o DECIMAL(5, 2),
    dg_kin_h2o DECIMAL(5, 2),
    ln_kf_den DECIMAL(5, 2),
    ln_ku_den DECIMAL(5, 2),
    mf DECIMAL(5, 2),
    mu DECIMAL(5, 2),
    cm DECIMAL(5, 2),
    authors VARCHAR(500),
    ref VARCHAR(500),
    publication_year INT,
    pmid VARCHAR(50),
    state VARCHAR(10),
    reversibility VARCHAR(10),
    review_date VARCHAR(50),
    additives VARCHAR(50),
    beta_t INT,
    secondary_structure VARCHAR(50),
    rsa DECIMAL(5, 2),
    dln_kf_h2o DECIMAL(5, 2),
    dln_ku_h2o DECIMAL(5, 2),
    dln_kf_den DECIMAL(5, 2),
    dln_ku_den DECIMAL(5, 2),
    ddgkin_h2o DECIMAL(5, 2),
    phif_h2o DECIMAL(5, 2),
    background_mutation_uniprot VARCHAR(50),
    background_mutation_pdb VARCHAR(50),
    reference_context VARCHAR(300),
    mutation_context VARCHAR(300),
    sequence_status VARCHAR(100)
);

CREATE TABLE IF NOT EXISTS protein(
    protein_id INT PRIMARY KEY AUTO_INCREMENT,
    protein_name VARCHAR(50),
    protein_source VARCHAR(50),
    protein_length DECIMAL(5, 2),
    uniprot_id VARCHAR(10),
	UNIQUE KEY uc_protein (protein_name, protein_source, protein_length, uniprot_id)
) ;


CREATE TABLE IF NOT EXISTS mutation(
    mutation_id INT PRIMARY KEY AUTO_INCREMENT,
    protein_id INT NOT NULL,
    mutated_chain VARCHAR(2),
    mutation_uniprot VARCHAR(20),
    mutation_pdb VARCHAR(20),
    UNIQUE KEY uc_mutation (protein_id, mutated_chain, mutation_uniprot, mutation_pdb),
    CONSTRAINT fk_protein_mutation
		FOREIGN KEY (protein_id)
        REFERENCES protein(protein_id)
);

CREATE TABLE IF NOT EXISTS structure(
    structure_id INT PRIMARY KEY AUTO_INCREMENT,
    protein_id INT NOT NULL,
    pdb_id VARCHAR(12),
    mutated_chain VARCHAR(2),
    secondary_structure VARCHAR(50),
    pfam VARCHAR(50),
    cath VARCHAR(50),
    ec_number VARCHAR(10),
    UNIQUE KEY uc_structure (protein_id, pdb_id, mutated_chain, secondary_structure, pfam, cath, ec_number), 
    CONSTRAINT fk_structure_mutation
		FOREIGN KEY (protein_id)
        REFERENCES protein(protein_id)
);

CREATE TABLE IF NOT EXISTS publication(
    publication_id INT PRIMARY KEY AUTO_INCREMENT,
    authors VARCHAR(500),
    ref VARCHAR(500),
    publication_year INT,
    pmid VARCHAR(50),
    review_date VARCHAR(50),
    CONSTRAINT uc_publication UNIQUE (publication_year, pmid, review_date)
);

CREATE TABLE IF NOT EXISTS experiment(
	experiment_id INT PRIMARY KEY AUTO_INCREMENT,
    protein_id INT NOT NULL,
    mutation_id INT NOT NULL,
    structure_id INT,
    publication_id INT NOT NULL,
    temperature DECIMAL(5, 2),
    ph DECIMAL(5, 2),
    buffer_name VARCHAR(50),
    buffer_conc VARCHAR(50),
    additives VARCHAR(50),
    kinetic_measurement VARCHAR(50),
    kinetic_method VARCHAR(50),
    state VARCHAR(10),
    reversibility VARCHAR(10),
    UNIQUE KEY uc_experiment(
		protein_id, 
        mutation_id, 
        structure_id, 
        publication_id, 
        temperature, 
        ph, 
        buffer_name,
        buffer_conc,
        additives,
        kinetic_measurement,
        kinetic_method,
        state,
        reversibility),
    CONSTRAINT fk_protein_experiment
		FOREIGN KEY (protein_id)
        REFERENCES protein(protein_id),
	CONSTRAINT fk_mutation_experiment
		FOREIGN KEY (mutation_id)
        REFERENCES mutation(mutation_id),
	CONSTRAINT fk_stucture_experiment
		FOREIGN KEY (structure_id)
        REFERENCES structure(structure_id),
	CONSTRAINT fk_publication_experiment
		FOREIGN KEY (publication_id)
        REFERENCES publication(publication_id)
);

CREATE TABLE IF NOT EXISTS measurement(
    measurement_id INT PRIMARY KEY AUTO_INCREMENT,
    experiment_id INT UNIQUE NOT NULL,
    ln_kf_h2o DECIMAL(5, 2),
    dln_kf_h2o DECIMAL(5, 2),
    ln_ku_h2o DECIMAL(5, 2),
    dln_ku_h2o DECIMAL(5, 2),
    dg_kin_h2o DECIMAL(5, 2),
    ddgkin_h2o DECIMAL(5, 2),
    mf DECIMAL(5, 2),
    mu DECIMAL(5, 2),
    beta_t DECIMAL(5, 2),
    CONSTRAINT fk_experiment_measurement
		FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id)
);

CREATE TABLE IF NOT EXISTS model_variable(
    variable_id INT PRIMARY KEY AUTO_INCREMENT,
    measurement_id INT NOT NULL UNIQUE,
    experiment_id INT NOT NULL,
    mutation_id INT NOT NULL,
    protein_id INT NOT NULL,
    reference_context VARCHAR(300),
    mutation_context VARCHAR(300),
    ln_kf_h2o DECIMAL(5, 2),
    dln_kf_h2o DECIMAL(5, 2),
    ln_ku_h2o DECIMAL(5, 2),
    dln_ku_h2o DECIMAL(5, 2),
    dg_kin_h2o DECIMAL(5, 2),
    ddgkin_h2o DECIMAL(5, 2),
    UNIQUE KEY uc_model(
		protein_id, 
        mutation_id, 
        measurement_id,
		ln_kf_h2o,
		dln_kf_h2o,
		ln_ku_h2o,
		dln_ku_h2o,
		dg_kin_h2o,
		ddgkin_h2o),
    CONSTRAINT fk_model_measurement
		FOREIGN KEY (measurement_id)
        REFERENCES measurement(measurement_id),
	CONSTRAINT fk_model_mutation
		FOREIGN KEY (mutation_id)
        REFERENCES mutation(mutation_id),
    CONSTRAINT fk_model_experiment
		FOREIGN KEY (experiment_id)
        REFERENCES experiment(experiment_id),
	CONSTRAINT fk_model_protein
		FOREIGN KEY (protein_id)
        REFERENCES protein(protein_id)
);


