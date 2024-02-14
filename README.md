# Comprehensive cheminformatics toolkit

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
    - `md_play_protocol.sh`: A shell script that automates the entire molecular dynamics process, from energy minimization to the production run.

4. **FreeEnergyLandscape**: Provides tools for analyzing and visualizing free energy landscapes from molecular dynamics simulations using collective variables (CVs). This directory includes a Python class designed to compute and plot free energy surfaces based on input data for two collective variables, such as angles and distances, facilitating the understanding of molecular mechanisms and energetics.
    - `FreeEnergyLandscape.py`: A comprehensive Python script that defines the `FreeEnergyLandscape` class. This class performs the loading of CV data, calculation of free energy using the Boltzmann inversion method, and visualization of free energy landscapes in 1D and 2D.
    - `GaussianKDE`: An optional, custom class implementation of the Gaussian Kernel Density Estimator (KDE) for estimating probability densities. This script is designed to replace dependencies on external libraries like SciPy for KDE calculations, making the toolkit more self-contained.


---

## Overview

Each script and input file in this repository is crafted to address specific needs in cheminformatics and molecular dynamics simulations. From data extraction and processing in ChEMBL to detailed molecular dynamics protocols and scaffold analysis, this toolkit serves as a comprehensive resource for researchers in the field.

The `FreeEnergyLandscape` directory introduces a set of tools essential for the analysis of energetic landscapes in molecular dynamics studies. By leveraging collective variables, these scripts offer insights into the thermodynamic stability and transition mechanisms of molecular systems, which are critical for understanding biological processes and material properties at the molecular level.


### How to use

Navigate to the desired directory and follow the instructions provided in each script or input file. Ensure you have the required dependencies installed, particularly for Python scripts.

To use the tools in the `FreeEnergyLandscape` directory, ensure you have a Python environment with NumPy and Matplotlib installed. Follow the instructions within each script for setting up your analysis. The `FreeEnergyLandscape.py` script can be directly executed with appropriate paths to your data files for CV1 and CV2, while `GaussianKDE` class is automatically utilized by the main script if present.

---

This repository is continuously updated with new tools and improved methodologies to keep up with the advancing field of cheminformatics and molecular dynamics simulations.
