USE kpro;

SELECT
    mv.variable_id,
    mv.measurement_id,
    e.experiment_id,

    -- Protein
    p.protein_id,
    p.uniprot_id,
    p.protein_name,
    p.protein_source,
    p.protein_length,

    -- Mutation
    mut.mutation_id,
    mut.mutated_chain,
    mut.mutation_uniprot,
    mut.mutation_pdb,

    -- Sequence-derived data
    mv.reference_context,
    mv.mutation_context,

    -- Target and related measurements
    m.ln_kf_h2o,
    m.dln_kf_h2o,
    m.ln_ku_h2o,
    m.dln_ku_h2o,
    m.dg_kin_h2o,
    m.ddgkin_h2o,
    m.mf,
    m.mu,
    m.beta_t,

    -- Experimental conditions
    e.temperature,
    e.ph,
    e.buffer_name,
    e.buffer_conc,
    e.additives,
    e.kinetic_measurement,
    e.kinetic_method,
    e.state,
    e.reversibility,

    -- Structure
    s.structure_id,
    s.pdb_id,
    s.secondary_structure,
    s.pfam,
    s.cath,
    s.ec_number,

    -- Publication provenance
    pub.publication_id,
    pub.ref,
    pub.pmid,
    pub.publication_year,
    pub.authors

FROM model_variable AS mv

JOIN measurement AS m
    ON m.measurement_id = mv.measurement_id

JOIN experiment AS e
    ON e.experiment_id = m.experiment_id
   AND e.experiment_id = mv.experiment_id

JOIN protein AS p
    ON p.protein_id = e.protein_id
   AND p.protein_id = mv.protein_id

JOIN mutation AS mut
    ON mut.mutation_id = e.mutation_id
   AND mut.mutation_id = mv.mutation_id

LEFT JOIN structure AS s
    ON s.structure_id = e.structure_id

JOIN publication AS pub
    ON pub.publication_id = e.publication_id

ORDER BY mv.variable_id;
