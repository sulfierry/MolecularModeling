# ChEMBL Database SQL Search and Post-Processing Scripts

This documentation describes two scripts utilized for querying the ChEMBL database using SQL and post-processing the obtained data.

## 1. SQL Script for ChEMBL Database Query
This script performs SQL queries in the ChEMBL database to retrieve kinase-related compounds with validated activity measurements.

### Script Overview:
- Connects to the ChEMBL database: `psql -U leon -d chembl_33`.
- Creates two tables:
  - `smile_kinase_manually_validated_kd_ki_ic50_10uM`: Contains compounds with IC50, Ki, and Kd values under 10µM and manually validated data.
  - `smile_kinase_all_compounds`: Includes all kinase-related compounds with activity measurements.
- Exports the results to TSV files.

### SQL:
```sql
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
```
