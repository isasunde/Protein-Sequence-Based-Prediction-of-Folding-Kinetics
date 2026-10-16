SELECT
	mv.variable_id,
	mv.measurement_id,
	e.experiment_id,

	p.protein_id,
	p.protein_name,
	p.uniprot_id,

	mut.mutation_id,
	mut.mutation_uniprot,
	mut.mutation_pdb,

	mv.reference_context,
	mv.mutation_context,

	m.ln_kf_h2o
    
FROM model_variable AS mv

JOIN measurement AS m
	ON m.measurement_id = mv.measurement_id

JOIN experiment AS e
	ON e.experiment_id = mv.experiment_id
	AND e.experiment_id = m.experiment_id

JOIN protein AS p
	ON p.protein_id = mv.protein_id
	AND p.protein_id = e.protein_id

JOIN mutation AS mut
	ON mut.mutation_id = mv.mutation_id
	AND mut.mutation_id = e.mutation_id

ORDER BY mv.measurement_id;