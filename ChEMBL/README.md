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

### SQL Commands:
```sql
-- Query for manually validated compounds with specific activity measurements
CREATE TABLE public.smile_kinase_manually_validated_kd_ki_ic50_10uM AS
SELECT DISTINCT [...]
\COPY public.smile_kinase_manually_validated_kd_ki_ic50_10uM TO '/home/leon/Desktop/ChEMBL_DATABSE/1_chembl_manually_validated/1_database/kinase_drug_info_all_manually_validated_IC50_Ki_kd_10uM.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t');

-- Query for all kinase-related compounds
CREATE TABLE public.smile_kinase_all_compounds AS
SELECT DISTINCT [...]
\COPY public.smile_kinase_all_compounds TO '/home/leon/Desktop/ChEMBL_DATABSE/1_chembl_manually_validated/1_database/kinase_all_compounds.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t');
```
