import pandas as pd

def process_kinase_group(group_name):
    special_cases = {
        "Enzyme": "Other", "Protein": "Other", "Other": "Other","Unclassified":"Other",
        "Atypical": "Other", "Kinase": "Other", "Nuclear e Other nuclear protein": "Nuclear protein",
        "Membrane": "Membrane receptor", "Voltage-gated potassium channel": "Potassium channel",
        "Tudor domain": "Tudor domain", "Bromodomain": "Bromodomain",
        "Transcription factor": "Transcription factor", "Plant homeodomain": "Plant homeodomain",
        "Protein Kinase": "Protein Kinase"
    }
    for key, value in special_cases.items():
        if group_name.startswith(key):
            return value
    return group_name.split()[0] if group_name else group_name

def filter_pkidb_chembl_ids(input_data, pkidb_file):
    # Load PKIDB data
    pkidb_data = pd.read_csv(pkidb_file, sep='\t')
    # Extract CHEMBL_IDs from PKIDB
    pkidb_ids = set(pkidb_data['CHEMBL_ID'])
    # Filter out entries in input_data that are in PKIDB
    filtered_data = input_data[~input_data['chembl_id'].isin(pkidb_ids)]
    return filtered_data

def abbreviate_kinase_groups(file_path, pkidb_file):
    # Carregar os dados do arquivo
    data = pd.read_csv(file_path, sep='\t')

    # Verificar se a coluna 'kinase_group' existe
    if 'kinase_group' in data.columns:
        # Aplicar a função process_kinase_group para processar os nomes dos grupos
        data['kinase_group'] = data['kinase_group'].apply(process_kinase_group)

    # Filter out CHEMBL_IDs from PKIDB
    data = filter_pkidb_chembl_ids(data, pkidb_file)

    return data

def main():
    file_path = './kinase_all_compounds.tsv'
    pkidb_file = './pkidb/pkidb_2024-03-18.tsv'
    processed_data = abbreviate_kinase_groups(file_path, pkidb_file)
    # Salvar os dados processados de volta para um arquivo ou visualizar como necessário
    processed_data.to_csv('./kinase_all_compounds_processed_kinase_groups.tsv', sep='\t', index=False)
    print("Processed data saved successfully.")

if __name__ == "__main__":
    main()
