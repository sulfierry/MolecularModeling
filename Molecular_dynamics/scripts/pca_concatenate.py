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

