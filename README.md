# Comprehensive Cheminformatics Toolkit

This repository is a diverse collection of tools and scripts aimed at facilitating various cheminformatics tasks. The directories are organized as follows:

1. **ChEMBL**: Contains SQL and Python scripts for querying and processing data from the ChEMBL database. It includes scripts for extracting kinase-targeted compounds and post-processing to remove salts and identify duplicates.
    - `descriptors.py`: Calculates molecular descriptors.
    - `cluster_by_similarity.py`: Clusters molecules based on structural similarity.
    - `histogram.py`: Generates histograms for molecular descriptors and similarity metrics.

2. **Scaffold_finder**: Focuses on scaffold analysis and similarity calculations in molecular datasets.
    - `substructure_scaffold_validator.py`: Validates substructures within molecules and highlights them.
    - `scaffold_finder.py`: Identifies molecular scaffolds and calculates similarity with a target structure.

3. **Molecular Dynamics**: Offers input files and protocols for molecular dynamics simulations, encompassing steps from energy minimization to long-term production runs.
    - Input parameters for different stages of molecular dynamics simulations.
    - Protocols detailing steps like pre-relaxation, relaxation, and production phases.

## Overview

Each script and input file in this repository is crafted to address specific needs in cheminformatics and molecular dynamics simulations. From data extraction and processing in ChEMBL to detailed molecular dynamics protocols and scaffold analysis, this toolkit serves as a comprehensive resource for researchers in the field.

### How to Use

Navigate to the desired directory and follow the instructions provided in each script or input file. Ensure you have the required dependencies installed, particularly for Python scripts.

---

This repository is continuously updated with new tools and improved methodologies to keep up with the advancing field of cheminformatics and molecular dynamics simulations.
