SELECT
	mv.variable_id,
	mv.measurement_id,
	mv.experiment_id,

	p.protein_id,
	p.protein_name,
	p.uniprot_id,

	mut.mutation_id,
	mut.mutation_uniprot,
	mut.mutation_pdb,

	mv.reference_context,
	mv.mutation_context,

	m.ln_kf_h2o,
	m.dln_kf_h2o,

	e.temperature,
	e.ph,
	e.buffer_name,
	e.buffer_conc,
	e.additives,
	e.kinetic_method,
	e.kinetic_measurement,
	e.state,
	e.reversibility,

	pub.publication_id,
	pub.pmid,
	pub.publication_year

FROM model_variable AS mv

JOIN measurement AS m
	ON m.measurement_id = mv.measurement_id
	AND m.experiment_id = mv.experiment_id

JOIN experiment AS e
	ON e.experiment_id = mv.experiment_id
	AND e.protein_id = mv.protein_id
	AND e.mutation_id = mv.mutation_id

JOIN protein AS p
	ON p.protein_id = mv.protein_id

JOIN mutation AS mut
	ON mut.mutation_id = mv.mutation_id

JOIN publication AS pub
	ON pub.publication_id = e.publication_id

ORDER BY mv.measurement_id;