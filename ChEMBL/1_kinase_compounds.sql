-- psql -U leon -d chembl_33
-- Criando a tabela kinase_group_activity_counts
CREATE TABLE kinase_group_activity_counts AS
SELECT kg.kinase_group, COUNT(*) AS count_kinase_group
FROM kinase_groups kg
JOIN assays a ON kg.tid = a.tid
JOIN activities act ON a.assay_id = act.assay_id
GROUP BY kg.kinase_group;


CREATE TABLE public.smile_kinase_all_compounds AS
SELECT DISTINCT
    d.chembl_id,
    cs.molregno,
    t.pref_name AS target_kinase,
    cs.canonical_smiles,
    act.standard_value,
    act.standard_type,
    act.pchembl_value,
    d.pref_name AS compound_name,
    t.organism AS organism
FROM
    compound_structures cs
JOIN
    activities act ON cs.molregno = act.molregno
JOIN
    assays a ON act.assay_id = a.assay_id
JOIN
    target_dictionary t ON a.tid = t.tid
LEFT JOIN
    molecule_dictionary d ON cs.molregno = d.molregno
WHERE
    t.pref_name LIKE '%kinase%' AND
    cs.canonical_smiles IS NOT NULL AND
    act.standard_type IN ('IC50', 'Ki', 'Kd') AND
    act.standard_value IS NOT NULL AND
    act.standard_units = 'nM' AND
    (act.data_validity_comment IS NULL OR act.data_validity_comment = 'Manually validated');

\COPY public.smile_kinase_all_compounds TO '/path/to/file/1_database/kinase_all_compounds.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t');
