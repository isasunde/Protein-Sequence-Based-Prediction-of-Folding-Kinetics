USE kpro;

INSERT IGNORE INTO protein (
    uniprot_id,
    protein_name,
    protein_source,
    protein_length
)
SELECT DISTINCT
    uniprot_id,
    protein_name,
    protein_source,
    protein_length
FROM placeholder
WHERE uniprot_id IS NOT NULL;

INSERT IGNORE INTO mutation (
    protein_id,
    mutated_chain,
    mutation_uniprot,
    mutation_pdb
)
SELECT DISTINCT
    p.protein_id,
    s.mutated_chain,
    s.mutation_uniprot,
    s.mutation_pdb
FROM placeholder AS s
JOIN protein AS p
    ON p.uniprot_id = s.uniprot_id
    AND p.protein_name = s.protein_name
    AND p.protein_source = s.protein_source
    AND p.protein_length = s.protein_length;
    
INSERT IGNORE INTO structure (
    protein_id,
    pdb_id,
    mutated_chain,
    secondary_structure,
    pfam,
    cath,
    ec_number
)
SELECT DISTINCT
    p.protein_id,
    s.pdb_id,
    s.mutated_chain,
    s.secondary_structure,
    s.pfam,
    s.cath,
    s.ec_number
FROM placeholder AS s
JOIN protein AS p
    ON p.uniprot_id = s.uniprot_id
    AND p.protein_name = s.protein_name
    AND p.protein_source = s.protein_source
    AND p.protein_length = s.protein_length;
    
INSERT IGNORE INTO publication (
    authors,
    ref,
    publication_year,
    pmid,
    review_date
)
SELECT DISTINCT
    authors,
    ref,
    publication_year,
    pmid,
    review_date
FROM placeholder;

INSERT IGNORE INTO experiment (
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
    reversibility
)
SELECT DISTINCT
    p.protein_id,
    m.mutation_id,
    s.structure_id,
    u.publication_id,
    l.temperature,
    l.ph,
    l.buffer_name,
    l.buffer_conc,
    l.additives,
    l.kinetic_measurement,
    l.kinetic_method,
    l.state,
    l.reversibility
FROM placeholder AS l
JOIN protein AS p
    ON p.uniprot_id = l.uniprot_id
    AND p.protein_name = l.protein_name
    AND p.protein_source = l.protein_source
    AND p.protein_length = l.protein_length
JOIN mutation AS m
	ON p.protein_id = m.protein_id
    AND l.mutated_chain = m.mutated_chain
    AND l.mutation_uniprot = m.mutation_uniprot
    AND l.mutation_pdb = m.mutation_pdb
JOIN structure AS s
	ON p.protein_id = s.protein_id
JOIN publication as u
	ON l.authors = u.authors
    AND l.ref = u.ref
    AND l.publication_year = u.publication_year
    AND l.pmid = u.pmid
    AND l.review_date = u.review_date;

INSERT IGNORE INTO measurement (
    experiment_id,
    ln_kf_h2o,
    dln_kf_h2o,
    ln_ku_h2o,
    dln_ku_h2o,
    dg_kin_h2o,
    ddgkin_h2o,
    mf,
    mu,
    beta_t
)
SELECT DISTINCT
    e.experiment_id,
    l.ln_kf_h2o,
    l.dln_kf_h2o,
    l.ln_ku_h2o,
    l.dln_ku_h2o,
    l.dg_kin_h2o,
    l.ddgkin_h2o,
    l.mf,
    l.mu,
    l.beta_t
FROM placeholder AS l
JOIN experiment AS e
    ON l.temperature = e.temperature
    AND l.ph = e.ph
    AND l.buffer_name = e.buffer_name
    AND l.buffer_conc = e.buffer_conc
    AND l.additives = e.additives
    AND l.kinetic_measurement = e.kinetic_measurement
    AND l.kinetic_method = e.kinetic_method
    AND l.state = e.state
    AND l.reversibility = e.reversibility;

TRUNCATE TABLE placeholder;