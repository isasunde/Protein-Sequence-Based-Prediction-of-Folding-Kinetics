USE kpro;

-- Row counts across the staging and normalized tables
SELECT 'placeholder' AS table_name, COUNT(*) AS row_count FROM placeholder
UNION ALL
SELECT 'protein', COUNT(*) FROM protein
UNION ALL
SELECT 'mutation', COUNT(*) FROM mutation
UNION ALL
SELECT 'structure', COUNT(*) FROM structure
UNION ALL
SELECT 'publication', COUNT(*) FROM publication
UNION ALL
SELECT 'experiment', COUNT(*) FROM experiment
UNION ALL
SELECT 'measurement', COUNT(*) FROM measurement
UNION ALL
SELECT 'model_variable', COUNT(*) FROM model_variable;


-- Experiments without measurements
SELECT e.experiment_id
FROM experiment AS e
LEFT JOIN measurement AS m
    ON m.experiment_id = e.experiment_id
WHERE m.measurement_id IS NULL;


-- Model variables whose measurement does not belong
-- to the same experiment
SELECT mv.variable_id
FROM model_variable AS mv
JOIN measurement AS m
    ON m.measurement_id = mv.measurement_id
WHERE mv.experiment_id <> m.experiment_id;


-- Experiments with invalid protein/mutation relationships
SELECT e.experiment_id
FROM experiment AS e
JOIN mutation AS m
    ON m.mutation_id = e.mutation_id
WHERE e.protein_id <> m.protein_id;


-- Experiments with a structure belonging to another protein
SELECT e.experiment_id
FROM experiment AS e
JOIN structure AS s
    ON s.structure_id = e.structure_id
WHERE e.protein_id <> s.protein_id;


SELECT
    l.record,
    l.uniprot_id,
    l.pdb_id,
    l.mutated_chain,
    l.pfam,
    l.cath,
    l.ec_number,
    COUNT(s.structure_id) AS matching_structures
FROM placeholder AS l
JOIN protein AS p
    ON p.protein_name <=> l.protein_name
   AND p.protein_source <=> l.protein_source
   AND p.protein_length <=> l.protein_length
   AND p.uniprot_id <=> l.uniprot_id
LEFT JOIN structure AS s
    ON s.protein_id = p.protein_id
   AND s.mutated_chain <=> l.mutated_chain
   AND s.pfam <=> l.pfam
   AND s.cath <=> l.cath
   AND s.ec_number <=> l.ec_number
GROUP BY
    l.record,
    l.uniprot_id,
    l.pdb_id,
    l.mutated_chain,
    l.pfam,
    l.cath,
    l.ec_number
HAVING COUNT(s.structure_id) > 1
ORDER BY matching_structures DESC;
