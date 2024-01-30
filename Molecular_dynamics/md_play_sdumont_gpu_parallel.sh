#!/bin/bash
#SBATCH --nodes=1                  # Número de nós
#SBATCH -p sequana_dockvs          # Fila (partition) a ser utilizada
#SBATCH -J molecular_dynamics      # Nome do job
#SBATCH --exclusive                # Utilização exclusiva dos nós
