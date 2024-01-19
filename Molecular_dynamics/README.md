# Molecular Dynamics Input Parameters

Este documento descreve os parâmetros comuns utilizados em simulações de dinâmica molecular, abrangendo as etapas de minimização, relaxamento e produção.

## Parâmetros Comuns em Todas as Etapas

- `imin`: Indica execução de minimização. `imin=1` para minimização de energia.
- `ntx`: Lê coordenadas do arquivo “inpcrd”. Opções suportadas: `1` (formato padrão sem velocidade inicial) e `2`.
- `irest`: Reinicia simulação. `irest=0` (padrão) não reinicia.
- `ntpr`: Imprime progresso a cada `ntpr` passos.
- `ntwx`: Grava coordenadas no arquivo `mdcrd` a cada `ntwx` passos. `ntwx=0` não grava.
- `cut`: Define corte não ligado em Angstroms. Valor comum: `8.0` para PME.

## Minimização
imin=1
ntx=1
irest=0
maxcyc=2000 # Máximo de ciclos de minimização
ncyc=1000 # Etapas de descida mais íngreme antes do gradiente conjugado
ntpr=100
ntwx=0
cut=10.0
