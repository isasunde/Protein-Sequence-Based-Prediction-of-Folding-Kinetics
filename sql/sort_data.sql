USE kpro;

START TRANSACTION;

# ============================================================
# 1. PROTEIN
# ============================================================
INSERT IGNORE INTO protein (
    protein_name,
    protein_source,
    protein_length,
    uniprot_id
)
SELECT DISTINCT
    protein_name,
    protein_source,
    protein_length,
    uniprot_id
FROM placeholder
WHERE uniprot_id IS NOT NULL;


# ============================================================
# 2. MUTATION
# ============================================================
INSERT IGNORE INTO mutation (
    protein_id,
    mutated_chain,
    mutation_uniprot,
    mutation_pdb
)
SELECT DISTINCT
    p.protein_id,
    l.mutated_chain,
    l.mutation_uniprot,
    l.mutation_pdb
FROM placeholder AS l
JOIN protein AS p
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id
WHERE l.uniprot_id IS NOT NULL;

# ============================================================
# 3. STRUCTURE
# ============================================================

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
    l.pdb_id,
    l.mutated_chain,
    l.secondary_structure,
    l.pfam,
    l.cath,
    l.ec_number
FROM placeholder AS l
JOIN protein AS p
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id
WHERE l.uniprot_id IS NOT NULL;

# ============================================================
# 4. PUBLICATION
# ============================================================
INSERT IGNORE INTO publication (
    authors,
    `ref`,
    publication_year,
    pmid,
    review_date
)
SELECT DISTINCT
    authors,
    `ref`,
    publication_year,
    pmid,
    review_date
FROM placeholder;

# ============================================================
# 5. EXPERIMENT
# ============================================================
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
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id

JOIN mutation AS m
    ON m.protein_id = p.protein_id
   AND m.mutated_chain <=> l.mutated_chain
   AND m.mutation_uniprot <=> l.mutation_uniprot
   AND m.mutation_pdb <=> l.mutation_pdb

JOIN structure AS s
    ON s.protein_id = p.protein_id
   AND s.pdb_id <=> l.pdb_id
   AND s.secondary_structure <=> l.secondary_structure
   AND s.mutated_chain <=> l.mutated_chain
   AND s.pfam <=> l.pfam
   AND s.cath <=> l.cath
   AND s.ec_number <=> l.ec_number

JOIN publication AS u
    ON u.publication_year <=> l.publication_year
   AND u.pmid <=> l.pmid
   AND u.review_date <=> l.review_date

WHERE l.uniprot_id IS NOT NULL;

# ============================================================
# 6. MEASUREMENT
# ============================================================

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

JOIN protein AS p
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id

JOIN mutation AS m
    ON m.protein_id = p.protein_id
   AND m.mutated_chain <=> l.mutated_chain
   AND m.mutation_uniprot <=> l.mutation_uniprot
   AND m.mutation_pdb <=> l.mutation_pdb

JOIN structure AS s
    ON s.protein_id = p.protein_id
   AND s.pdb_id <=> l.pdb_id
   AND s.mutated_chain <=> l.mutated_chain
   AND s.secondary_structure <=> l.secondary_structure
   AND s.pfam <=> l.pfam
   AND s.cath <=> l.cath
   AND s.ec_number <=> l.ec_number

JOIN publication AS u
    ON u.publication_year <=> l.publication_year
   AND u.pmid <=> l.pmid
   AND u.review_date <=> l.review_date

JOIN experiment AS e
    ON e.protein_id = p.protein_id
   AND e.mutation_id = m.mutation_id
   AND e.structure_id = s.structure_id
   AND e.publication_id = u.publication_id
   AND e.temperature <=> l.temperature
   AND e.ph <=> l.ph
   AND e.buffer_name <=> l.buffer_name
   AND e.buffer_conc <=> l.buffer_conc
   AND e.kinetic_measurement <=> l.kinetic_measurement
   AND e.kinetic_method <=> l.kinetic_method;

# ============================================================
# 7. MODEL VARIABLE
# ============================================================

INSERT IGNORE INTO model_variable (
    measurement_id,
    experiment_id,
    mutation_id,
    protein_id,
    reference_context,
    mutation_context,
    ln_kf_h2o,
    dln_kf_h2o,
    ln_ku_h2o,
    dln_ku_h2o,
    dg_kin_h2o,
    ddgkin_h2o
)
SELECT DISTINCT
    m.measurement_id,
    e.experiment_id,
    e.mutation_id,
    e.protein_id,
    l.reference_context,
    l.mutation_context,
    m.ln_kf_h2o,
    m.dln_kf_h2o,
    m.ln_ku_h2o,
    m.dln_ku_h2o,
    m.dg_kin_h2o,
    m.ddgkin_h2o
FROM placeholder AS l

JOIN protein AS p
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id

JOIN mutation AS u
    ON u.protein_id = p.protein_id
   AND u.mutated_chain <=> l.mutated_chain
   AND u.mutation_uniprot <=> l.mutation_uniprot
   AND u.mutation_pdb <=> l.mutation_pdb

JOIN structure AS s
    ON s.protein_id = p.protein_id
   AND s.pdb_id <=> l.pdb_id
   AND s.secondary_structure <=> l.secondary_structure
   AND s.mutated_chain <=> l.mutated_chain
   AND s.pfam <=> l.pfam
   AND s.cath <=> l.cath
   AND s.ec_number <=> l.ec_number

JOIN publication AS pub
    ON pub.publication_year <=> l.publication_year
   AND pub.pmid <=> l.pmid
   AND pub.review_date <=> l.review_date

JOIN experiment AS e
    ON e.protein_id = p.protein_id
   AND e.mutation_id = u.mutation_id
   AND e.structure_id = s.structure_id
   AND e.publication_id = pub.publication_id
   AND e.temperature <=> l.temperature
   AND e.ph <=> l.ph
   AND e.buffer_name <=> l.buffer_name
   AND e.buffer_conc <=> l.buffer_conc
   AND e.kinetic_measurement <=> l.kinetic_measurement
   AND e.kinetic_method <=> l.kinetic_method

JOIN measurement AS m
    ON m.experiment_id = e.experiment_id

WHERE l.uniprot_id IS NOT NULL;

COMMIT;