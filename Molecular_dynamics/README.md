# Molecular Dynamics Play Protocol (`md_play_protocol.sh`)

This shell script automates the execution of various stages in molecular dynamics simulations using AMBER software. It orchestrates the transition from energy minimization to production runs, handling each step with appropriate input parameters and file management.

## Overview

The script is designed to sequentially process multiple stages of molecular dynamics:

1. **Energy Minimization**: Reduces the energy of the system to a local minimum.
2. **Relaxation**: Involves several phases to gradually relax the system, each with specific constraints and timeframes.
3. **Pre-Production**: Equilibrates the system under defined conditions before the main production run.
4. **Production**: The final stage where the actual molecular dynamics simulation is carried out.

## How it Works

- The script defines arrays for input files, restart files, trajectories, and outputs corresponding to each stage.
- It uses `pmemd.cuda` for GPU-accelerated processing and `sander.MPI` for parallel execution.
- The script iterates through each stage, executing the AMBER tools with the specified parameters.
- Conditional logic is applied to handle different requirements at each simulation stage.

### Key Stages and Parameters

- **Minimization (1000 steps)**: Uses the `0_minimization.in` input file.
- **Relaxation Phases (total 1 ns)**:
  - **Phase 1**: 300 ps pre-relaxation NPT with constraints on protein and ligand.
  - **Phase 2**: 300 ps relaxation NPT with protein constraints.
  - **Phase 3**: 200 ps relaxation NPT, no constraints on side chains within 5 Å of the ligand.
  - **Phase 4**: 200 ps relaxation NPT, no constraints on residues within 5 Å of the ligand.
- **Pre-Production (5 ns)**: Prepares the system for the production run.
- **Production (100 ns)**: The main simulation phase, capturing the dynamics of the molecular system.

## Usage

Run the script in a shell environment with access to AMBER tools. Ensure the input files and `prmtop` file are correctly set up.

```bash
./md_play_protocol.sh
```

---

This comprehensive simulation protocol ensures that the system is adequately prepared and equilibrated at each stage before proceeding to the long-term dynamic analysis in the production phase.


