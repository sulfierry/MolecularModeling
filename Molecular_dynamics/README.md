# Molecular dynamics folder overview

This directory is dedicated to the setup and execution of molecular dynamics simulations. It contains organized subfolders and scripts that streamline the parametrization and simulation processes. Below is an overview of each component within this folder:

## Subfolders

1. **`ligand_parametrization`**:
   - Contains scripts and input files related to the parametrization of ligands. This process is essential for preparing ligands for accurate simulation within the molecular dynamics framework.

2. **`parametrization (protein)`**:
   - This folder includes scripts and tools for the parametrization of the entire molecular system, including proteins, ligands, and solvent models. It typically involves setting up force fields, system solvation, and ion placements.
    
3. **`molecular_dynamics_input`**:
   - Houses input files such as initial coordinates and configurations required for running molecular dynamics simulations. These inputs are foundational for initiating and accurately running simulations.


## Shell script

- **`md_play_protocoll.sh`**:
   - A shell script designed to automate the process of running a series of molecular dynamics simulations. The script organizes the workflow into distinct stages, such as minimization, relaxation, pre-production, and production phases. It ensures that each step is executed sequentially and efficiently, utilizing the input files and parameters defined in the above folders.

### Script functions
   - Each stage of the simulation is defined with specific input parameters and configurations, ensuring a comprehensive exploration of the molecular system's dynamics.
   - The script supports both GPU-accelerated (`pmemd.cuda`) and parallel CPU (`mpirun -np 48 sander.MPI`) executions, making it adaptable to various computational resources.
   - Output files from each stage are systematically saved, providing a detailed trajectory of the simulation process.

#### Overview

The script is designed to sequentially process multiple stages of molecular dynamics:

1. **Energy Minimization**: Reduces the energy of the system to a local minimum.
2. **Relaxation**: Involves several phases to gradually relax the system, each with specific constraints and timeframes.
3. **Pre-Production**: Equilibrates the system under defined conditions before the main production run.
4. **Production**: The final stage where the actual molecular dynamics simulation is carried out.

#### How it works

- The script defines arrays for input files, restart files, trajectories, and outputs corresponding to each stage.
- It uses `pmemd.cuda` for GPU-accelerated processing and `sander.MPI` for parallel execution.
- The script iterates through each stage, executing the AMBER tools with the specified parameters.
- Conditional logic is applied to handle different requirements at each simulation stage.

##### Key stages and parameters

- **Minimization (1000 steps)**: Uses the `0_minimization.in` input file.
- **Relaxation Phases (total 1 ns)**:
  - **Phase 1**: 300 ps pre-relaxation NPT with constraints on protein and ligand.
  - **Phase 2**: 300 ps relaxation NPT with protein constraints.
  - **Phase 3**: 200 ps relaxation NPT, no constraints on side chains within 5 Å of the ligand.
  - **Phase 4**: 200 ps relaxation NPT, no constraints on residues within 5 Å of the ligand.
- **Pre-Production (5 ns)**: Prepares the system for the production run.
- **Production (100 ns)**: The main simulation phase, capturing the dynamics of the molecular system.

## Usage

- To run the molecular dynamics protocol, navigate to this directory and execute the `md_play_protocoll.sh` script. Ensure that all prerequisite files and parameters are correctly set up in the respective subfolders.

Run the script in a shell environment with access to AMBER tools. Ensure the input files and `prmtop` file are correctly set up.

```bash
./md_play_protocol.sh
```
- This comprehensive simulation protocol ensures that the system is adequately prepared and equilibrated at each stage before proceeding to the long-term dynamic analysis in the production phase.


## Utilização Paralela das GPUs

O script `md_play_sdumont_gpu_parallel.sh` foi desenvolvido para otimizar a utilização dos recursos computacionais no SDumont, permitindo a execução paralela dos jobs em todas as GPUs disponíveis no atual nó do cluster. Isso é particularmente útil para simulações de dinâmica molecular que podem se beneficiar significativamente do paralelismo em GPU.

### Alocação de GPUs

Para gerenciar qual GPU é usada por cada tarefa, o script define a variável de ambiente `CUDA_VISIBLE_DEVICES` usando a seguinte lógica:

```bash
export CUDA_VISIBLE_DEVICES=$(( (REPLICA - 1) % 4 ))

# Esta expressão calcula um índice de 0 a 3, correspondente a cada uma das 4 GPUs,
# baseando-se no número sequencial da réplica. Dessa forma, cada réplica é alocada a uma
# GPU específica, permitindo a distribuição equitativa das tarefas pelas GPUs disponíveis.

```

### Execução Paralela
O script emprega um loop que submete conjuntos de até 4 réplicas em paralelo, cada uma utilizando uma GPU diferente. Após a submissão de um conjunto de réplicas, o script aguarda a conclusão dessas tarefas antes de proceder com o próximo conjunto. Isso assegura que cada tarefa seja executada na sua respectiva GPU sem sobreposição ou conflito de recursos.

```bash
# Número total de réplicas
total_replicas=10
# Número de GPUs disponíveis
num_gpus=4

...

# Loop para executar as réplicas em grupos de 4
for (( i=1; i<=total_replicas; i+=num_gpus )); do
    for (( j=i; j<i+num_gpus && j<=total_replicas; j++ )); do
        run_replica $j &
    done
    wait
done
```

### Vantagens
A execução paralela de tarefas utilizando o script run.sh oferece várias vantagens:

- **Eficiência**: Maximiza a utilização dos recursos de GPU, reduzindo o tempo total de processamento.
- **Escalabilidade**: Facilmente ajustável para diferentes quantidades de réplicas ou configurações de sistema.
- **Simplicidade**: Permite aos usuários submeterem múltiplas réplicas sem necessidade de gerenciamento manual dos recursos de GPU.

---

This organized structure within the "Molecular Dynamics" folder facilitates a streamlined and efficient approach to conducting molecular dynamics simulations, catering to both novice and experienced users in the field of computational chemistry and molecular modeling.



