# Substructure and Scaffold Validator

The script `substructure_scaffold_validator.py` focuses on validating the presence of specific substructures within larger molecular structures. It utilizes RDKit for molecular manipulation and visualization.

Key features:
- Converts SMILES strings to RDKit molecule objects.
- Defines a target substructure in SMILES format.
- Checks if the substructure exists within the full molecule.
- Highlights the substructure within the full molecular structure and saves the image.
- Displays both the complete molecule and the substructure for visual validation.

![Alt text da image](https://github.com/sulfierry/MolecularModelingTools/blob/main/Scaffold_finder/highlighted_substructure.png)

This script is particularly useful for confirming the presence of key functional groups or motifs in larger molecules.


# Scaffold Finder and Similarity Calculator

`scaffold_finder.py` is designed to identify molecular scaffolds and calculate similarity with a target structure. It utilizes RDKit for molecular fingerprinting and similarity calculations.

Primary functionalities:
- Extracts Murcko scaffolds from molecules based on their SMILES strings.
- Calculates Tanimoto similarity between each molecule and a predefined target SMILES.
- Processes a list of molecules in parallel to efficiently handle large datasets.
- Adds scaffold information and similarity scores to the existing dataset.
- Filters and saves the dataset based on similarity thresholds.



This script aids in the identification of common scaffolds and the assessment of molecular similarity, which is crucial in cheminformatics and drug discovery processes.
