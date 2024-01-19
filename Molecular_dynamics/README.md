# Molecular Dynamics Folder Overview

This directory is dedicated to the setup and execution of molecular dynamics simulations. It contains organized subfolders and scripts that streamline the parametrization and simulation processes. Below is an overview of each component within this folder:

## Subfolders

1. **`ligand_parametrization`**:
   - Contains scripts and input files related to the parametrization of ligands. This process is essential for preparing ligands for accurate simulation within the molecular dynamics framework.

2. **`molecular_dynamics_input`**:
   - Houses input files such as initial coordinates and configurations required for running molecular dynamics simulations. These inputs are foundational for initiating and accurately running simulations.

3. **`parametrization`**:
   - This folder includes scripts and tools for the parametrization of the entire molecular system, including proteins, ligands, and solvent models. It typically involves setting up force fields, system solvation, and ion placements.

## Shell Script

- **`md_play_protocoll.sh`**:
   - A shell script designed to automate the process of running a series of molecular dynamics simulations. The script organizes the workflow into distinct stages, such as minimization, relaxation, pre-production, and production phases. It ensures that each step is executed sequentially and efficiently, utilizing the input files and parameters defined in the above folders.

### Script Functions
   - Each stage of the simulation is defined with specific input parameters and configurations, ensuring a comprehensive exploration of the molecular system's dynamics.
   - The script supports both GPU-accelerated (`pmemd.cuda`) and parallel CPU (`mpirun -np 48 sander.MPI`) executions, making it adaptable to various computational resources.
   - Output files from each stage are systematically saved, providing a detailed trajectory of the simulation process.

## Usage

- To run the molecular dynamics protocol, navigate to this directory and execute the `md_play_protocoll.sh` script. Ensure that all prerequisite files and parameters are correctly set up in the respective subfolders.

This organized structure within the "Molecular Dynamics" folder facilitates a streamlined and efficient approach to conducting molecular dynamics simulations, catering to both novice and experienced users in the field of computational chemistry and molecular modeling.


## Molecular Dynamics Play Protocol (`md_play_protocol.sh`)

This shell script automates the execution of various stages in molecular dynamics simulations using AMBER software. It orchestrates the transition from energy minimization to production runs, handling each step with appropriate input parameters and file management.

### Overview

The script is designed to sequentially process multiple stages of molecular dynamics:

1. **Energy Minimization**: Reduces the energy of the system to a local minimum.
2. **Relaxation**: Involves several phases to gradually relax the system, each with specific constraints and timeframes.
3. **Pre-Production**: Equilibrates the system under defined conditions before the main production run.
4. **Production**: The final stage where the actual molecular dynamics simulation is carried out.

### How it Works

- The script defines arrays for input files, restart files, trajectories, and outputs corresponding to each stage.
- It uses `pmemd.cuda` for GPU-accelerated processing and `sander.MPI` for parallel execution.
- The script iterates through each stage, executing the AMBER tools with the specified parameters.
- Conditional logic is applied to handle different requirements at each simulation stage.

#### Key Stages and Parameters

- **Minimization (1000 steps)**: Uses the `0_minimization.in` input file.
- **Relaxation Phases (total 1 ns)**:
  - **Phase 1**: 300 ps pre-relaxation NPT with constraints on protein and ligand.
  - **Phase 2**: 300 ps relaxation NPT with protein constraints.
  - **Phase 3**: 200 ps relaxation NPT, no constraints on side chains within 5 Å of the ligand.
  - **Phase 4**: 200 ps relaxation NPT, no constraints on residues within 5 Å of the ligand.
- **Pre-Production (5 ns)**: Prepares the system for the production run.
- **Production (100 ns)**: The main simulation phase, capturing the dynamics of the molecular system.

### Usage

Run the script in a shell environment with access to AMBER tools. Ensure the input files and `prmtop` file are correctly set up.

```bash
./md_play_protocol.sh
```

---

This comprehensive simulation protocol ensures that the system is adequately prepared and equilibrated at each stage before proceeding to the long-term dynamic analysis in the production phase.


