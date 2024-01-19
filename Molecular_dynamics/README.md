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
```
imin=1
ntx=1
irest=0
maxcyc=2000 # Máximo de ciclos de minimização
ncyc=1000 # Etapas de descida mais íngreme antes do gradiente conjugado
ntpr=100
ntwx=0
cut=10.0
```

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
```
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
```

### Parâmetros Explicados

- `dt`: Passo de tempo. Máximo recomendado: `0.002` com SHAKE.
- `ntf` e `ntc`: Avaliação de força e restrições de comprimento de ligação com SHAKE.
- `tempi`: Temperatura inicial.
- `temp0`: Temperatura de referência. Padrão: `300K`.
- `ntwx` e `ntwr`: Frequência de gravação de coordenadas e arquivo de reinício.
- `ntb`: Controle de limites periódicos. Padrão para pressão constante: `ntb=2`.
- `ntp`: Dinâmica de pressão constante. `ntp=1` para escala de posição isotrópica.
- `ntt`: Escalonamento de temperatura. `ntt=2` para conjunto canônico.  Observe que a configuração `ntt=0` corresponde ao ensemble microcanônico `(NVE)` (que deve se aproximar do canônico para grandes números de graus de liberdade). Alguns aspectos do "conjunto de acoplamento fraco" `(ntt=1)` foram examinados e aproximadamente interpolados entre os conjuntos microcanônico e canônico. As opções `ntt=2` e `3` correspondem ao conjunto canônico (T constante).
- `ntt=2` pois definimos o esquema de acoplamento de temperatura do tipo `Andersen`, no qual "colisões" imaginárias randomizam as velocidades para uma distribuição correspondente a `temp0` a cada `vrand` passos. Observe que entre essas "colisões massivas", a dinâmica é newtoniana. Portanto, as funções de correlação de tempo (etc.) podem ser calculadas nessas seções e os resultados calculados em uma distribuição canônica inicial. Observe também que uma taxa de colisão muito alta (um valor muito pequeno de `vrand`) diminuirá a velocidade com que as moléculas exploram o espaço de configuração, enquanto uma taxa muito baixa significa que a distribuição canônica de energias será amostrada lentamente. Uma discussão sobre essa taxa é fornecida por `Andersen`. Observe que esta opção não é equivalente ao termostato original descrito por `Andersen`.

- `nmropt`: Controle de parâmetros NMR e variações. `nmropt=1` para restrições e alterações.

- `ig`: Semente para gerador de números pseudo-aleatórios. `-1` para semente baseada em data/hora.
- `iwrap`: Empacotamento de coordenadas em caixa primária. `iwrap=1` para empacotamento. Isso significa que, para cada molécula, sua imagem periódica mais próxima do meio da "caixa primária" (com coordenadas x entre 0 e a, coordenadas y entre 0 e b e coordenadas z entre 0 e c) será aquela escrita para o arquivo de saída. Isso geralmente faz com que as estruturas resultantes pareçam visualmente melhores, mas não tem efeito sobre a energia ou as forças. Executar tal empacotamento, no entanto, pode atrapalhar a difusão e outros cálculos.
- Para execuções muito longas, pode ser necessário definir `iwrap = 1` para evitar que a saída de coordenadas transborde a trajetória e reinicie os formatos de arquivo, especialmente se as trajetórias forem escritas no formato `ASCII` em vez de `NetCDF` (consulte também a opção `ioutfm`).
- `ntr`: Restrição de átomos no espaço cartesiano utilizando potencial harmônico. `ntr=1` com `restraintmask` definido.

---

Esse formato proporciona uma estrutura clara e concisa, facilitando a compreensão dos parâmetros e suas funções nas simulações de dinâmica molecular.
