-- psql -U leon -d chembl_33

CREATE TABLE public.smile_kinase_manually_validated_kd_ki_ic50_10uM AS
SELECT DISTINCT
    d.chembl_id,
    cs.molregno,
    t.pref_name AS target_kinase,
    cs.canonical_smiles,
    act.standard_value,
    act.standard_type,
    act.pchembl_value,
    d.pref_name AS compound_name,
    o.organism AS organism
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
LEFT JOIN
    component_sequences o ON t.component_id = o.component_id
WHERE
    t.pref_name LIKE '%kinase%' AND
    cs.canonical_smiles IS NOT NULL AND
    act.standard_type IN ('IC50', 'Ki', 'Kd') AND
    act.standard_value IS NOT NULL AND
    act.standard_units = 'nM' AND
    act.standard_value < 10000 AND
    act.standard_relation = '=' AND
    (act.data_validity_comment IS NULL OR act.data_validity_comment = 'Manually validated');

\COPY public.smile_kinase_manually_validated_kd_ki_ic50_10uM TO '/home/leon/Desktop/ChEMBL_DATABSE/1_chembl_manually_validated/1_database/kinase_drug_info_all_manually_validated_IC50_Ki_kd_10uM.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t');


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

\COPY public.smile_kinase_all_compounds TO '/home/leon/Desktop/ChEMBL_DATABSE/1_chembl_manually_validated/1_database/kinase_all_compounds.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t');
