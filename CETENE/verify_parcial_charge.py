def soma_cargas_mol2(arquivo_mol2):
    carga_total = 0.0
    with open(arquivo_mol2, 'r') as file:
        for line in file:
            if line.startswith("@<TRIPOS>ATOM"):
                break
        for line in file:
            if line.startswith("@<TRIPOS>"):
                break
            partes = line.split()
            if partes:
                carga = float(partes[-1])
                carga_total += carga
    
    return carga_total
