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


> `maxcyc`: Máximo de ciclos de minimização. `maxcyc=2000` para minimizador `xmin` com método TNCG.

### LMOD Procedure

- Utiliza autovetores de modos vibracionais de baixa frequência.
- Inicia com modelo molecular minimizado energeticamente.
- Usa ARPACK para modos vibracionais.
- `Hv` calculado via diferenças finitas: `Hv = [∇(xmin+h)−∇(xmin)]/h`.
- Reutilização de autovetores para eficiência em LMOD.

### XMIN
float xmin( float func(), int natm, float x[], float g[], float ene, float grms_out, struct xmod_opt xo);


- `xmin()`: Retorna energia minimizada e atualiza coordenadas para a conformação de energia mínima.

## Configurações para Relaxamento, Pré-produção e Produção

imin=0
dt=0.002
ntf=2
ntc=2
temp0=300.0
ntpr=100
ntwx=2000
ntwr=2000
cut=10.0
ntb=2
ntp=1
ntt=2
nmropt=1
ig=-1
iwrap=1
ntr=1
restraintmask="X"


### Parâmetros Explicados

- `dt`: Passo de tempo. Máximo recomendado: `0.002` com SHAKE.
- `ntf` e `ntc`: Avaliação de força e restrições de comprimento de ligação com SHAKE.
- `temp0`: Temperatura de referência. Padrão: `300K`.
- `ntwx` e `ntwr`: Frequência de gravação de coordenadas e arquivo de reinício.
- `ntb`: Controle de limites periódicos. Padrão para pressão constante: `ntb=2`.
- `ntp`: Dinâmica de pressão constante. `ntp=1` para escala de posição isotrópica.
- `ntt`: Escalonamento de temperatura. `ntt=2` para conjunto canônico.
- `nmropt`: Controle de parâmetros NMR e variações. `nmropt=1` para restrições e alterações.
- `ig`: Semente para gerador de números pseudo-aleatórios. `-1` para semente baseada em data/hora.
- `iwrap`: Empacotamento de coordenadas em caixa primária. `iwrap=1` para empacotamento.
- `ntr`: Restrição de átomos no espaço cartesiano. `ntr=1` com `restraintmask` definido.

---

Esse formato proporciona uma estrutura clara e concisa, facilitando a compreensão dos parâmetros e suas funções nas simulações de dinâmica molecular.
