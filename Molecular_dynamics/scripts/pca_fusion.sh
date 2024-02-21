#!/bin/bash

# apos remover as aguas e processar atraves da classe "MDTraj"
# unir os valores de pca em apenas um arquivo (pca1_all.tsv e pca2_all.tsv)

# Definindo o diretório onde as pastas 1 a 10 estão localizadas
# Substitua "/caminho/para/o/diretorio" pelo caminho real em seu sistema
base_dir="./"

# Criando a pasta cv_all_together se ela não existir
mkdir -p "$base_dir/cv_all_together"

# Inicializando arquivos finais
> "$base_dir/cv_all_together/cv_pca1_all.tsv"
> "$base_dir/cv_all_together/cv_pca2_all.tsv"

# Loop através das pastas de 1 a 10
for i in {1..10}; do
    # Concatenando o conteúdo de cv_pca1.tsv para cv_pca1_all.tsv
    if [ $i -eq 1 ]; then
        cat "$base_dir/$i/cv_pca1.tsv" >> "$base_dir/cv_all_together/cv_pca1_all.tsv"
    else
        tail -n +2 "$base_dir/$i/cv_pca1.tsv" >> "$base_dir/cv_all_together/cv_pca1_all.tsv"
    fi
    # Concatenando o conteúdo de cv_pca2.tsv para cv_pca2_all.tsv
    if [ $i -eq 1 ]; then
        cat "$base_dir/$i/cv_pca2.tsv" >> "$base_dir/cv_all_together/cv_pca2_all.tsv"
    else
        tail -n +2 "$base_dir/$i/cv_pca2.tsv" >> "$base_dir/cv_all_together/cv_pca2_all.tsv"
    fi
done
