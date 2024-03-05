import pandas as pd
from tqdm.auto import tqdm
import concurrent.futures
import os


# Função para ler e processar um arquivo em chunks
def process_file_chunkwise(filename, start_idx):
    chunk_size = 10240  # Defina o tamanho do chunk aqui
    chunks = []
    for chunk in pd.read_csv(filename, sep='\t', chunksize=chunk_size):
        chunk.index = range(start_idx, start_idx + len(chunk))  # Reindexa
        chunks.append(chunk)
        start_idx += len(chunk)
    return pd.concat(chunks)

# Função para processar arquivos em paralelo e concatená-los
def concatenate_files(file_list):
    # Determina os índices de início para cada arquivo
    start_indices = [0]
    for filename in file_list[:-1]:
        start_indices.append(start_indices[-1] + sum(1 for line in open(filename)) - 1)
    
    # Processa os arquivos em paralelo
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_file_chunkwise, file, start_idx) for file, start_idx in zip(file_list, start_indices)]
        results = [f.result() for f in tqdm(concurrent.futures.as_completed(futures), total=len(file_list))]
    
    # Concatena os resultados
    final_df = pd.concat(results)
    return final_df

