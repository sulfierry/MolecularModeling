# ChEMBL Database SQL Search and Post-Processing Scripts

This documentation describes two scripts utilized for querying the ChEMBL database using SQL and post-processing the obtained data.

## 1. SQL Script for ChEMBL Database Query - kinase_compounds.sql
This script performs SQL queries in the ChEMBL database to retrieve kinase-related compounds with validated activity measurements.

### Script Overview:
- Connects to the ChEMBL database: `psql -U leon -d chembl_33`.
- Creates two tables:
  - `smile_kinase_manually_validated_kd_ki_ic50_10uM`: Contains compounds with IC50, Ki, and Kd values under 10µM and manually validated data.
  - `smile_kinase_all_compounds`: Includes all kinase-related compounds with activity measurements.
- Exports the results to TSV files.

#### Script Execution
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

## 2. Python Post-Processing Script - remove_redundance.py:

The `verify_partial_charge.py` Python script is employed for post-processing kinase-related compound data, focusing on removing redundancies and identifying salt-free compounds.

#### Key Functions of the Script:

1. **`remove_salts_and_identify(smiles)`**: 
   - Removes salts from the SMILES of molecules and identifies if the original molecule contained salts.
   - Returns the salt-free SMILES and a boolean indicating the presence of salt.

2. **`redundance_remov(file_path)`**: 
   - Loads data from the TSV file and identifies molecules containing salts.
   - Groups by canonical SMILES to get lists of target kinases.
   - Selects the row with the lowest activity value for each SMILES, removing redundancies.
   - Adds additional columns for detailed analyses.

3. **`create_salt_free_output(processed_data)`**: 
   - Creates a salt-free version of the processed data.

4. **`save_positive_negative_files(salt_free_data, output_directory)`**: 
   - Saves two separate TSV files: one for compounds with positive activity and another for negative, based on the activity value.

#### Script Execution

```python
def main():
    input_file_path = './kinase_all_compounds_updated.tsv'
    output_file_path = './nr_kinase_all_compounds_updated.tsv'
    salt_free_output_path = './nr_kinase_all_compounds_salt_free.tsv'
    output_directory = '.'  # Directory for 'positive.tsv' and 'negative.tsv'

    # Process and remove redundancies
    processed_data = redundance_remov(input_file_path)
    processed_data.to_csv(output_file_path, sep='\t', index=False, na_rep='')

    # Create salt-free output
    salt_free_data = create_salt_free_output(processed_data)
    salt_free_data.to_csv(salt_free_output_path, sep='\t', index=False, na_rep='')

    # Save 'positive.tsv' and 'negative.tsv' files using salt-free data
    save_positive_negative_files(salt_free_data, output_directory)

    print(f"Processed file saved at: {output_file_path}")
    print(f"Salt-free processed file saved at: {salt_free_output_path}")
    print("Files 'positive.tsv' and 'negative.tsv' have been created.")

if __name__ == "__main__":
    main()
```

# Molecular Descriptors Calculation

This script, `descriptors.py`, is designed to calculate various molecular descriptors from SMILES strings. These descriptors include molecular weight (MW), LogP, Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), Topological Polar Surface Area (TPSA), and Number of Rotatable Bonds (NRB). The script reads the SMILES strings from a given file, computes these descriptors using RDKit, and saves the results.

Key functionalities include:
- Reading SMILES strings from a TSV file.
- Calculating molecular descriptors for each molecule.
- Saving the computed descriptors to a TSV file.
- Generating histograms to visualize the distribution of these descriptors.

![Alt text da image](https://github.com/sulfierry/MolecularModelingTools/blob/main/ChEMBL/nr_chembl_pkidb_descriptors.png)


This script is an essential tool in the pre-processing of molecular data for cheminformatics analyses.


# Molecular Clustering Based on Similarity

`cluster_by_similarity.py` focuses on grouping molecules based on structural similarities. It uses the SMILES representation of molecules, converting these structures into molecular fingerprints. The script then performs clustering based on a predefined similarity threshold.

Highlights of this script include:
- Conversion of SMILES to molecular fingerprints using RDKit.
- All-against-all comparison of fingerprints for similarity calculation.
- Clustering of molecules based on Tanimoto similarity threshold.
- Visualization of clusters in 2D space using t-SNE.

![Alt text da image](https://github.com/sulfierry/MolecularModelingTools/blob/main/ChEMBL/tsne_similarity.png)

This script offers an insightful approach to understand molecular diversity and similarity in large datasets.

# Histograms for Molecular Similarity and Distance

The `histogram.py` script compares molecular fingerprints in an all-against-all approach, computing similarity and distance measures for every pair of molecules in the dataset. It processes fingerprints in batches to efficiently manage memory and computational resources.

Key features of the script:
- Loading molecular fingerprints from SMILES strings.
- Computing similarity using Tanimoto, Dice, and Cosine metrics.
- Computing distance using Hamming and Manhattan metrics.
- Generating histograms to display the distributions of these similarity and distance measures.

This script is particularly useful for assessing the overall structural diversity and similarity within a molecular dataset.

