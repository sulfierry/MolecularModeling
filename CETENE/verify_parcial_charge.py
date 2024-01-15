import matplotlib.pyplot as plt

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

# Substitua 'molecula.mol2' pelo caminho do seu arquivo .mol2
carga_total = soma_cargas_mol2('output_test.mol2')

if carga_total > 0:
    print(f"A carga total é positiva: {carga_total:.4f}")
elif carga_total < 0:
    print(f"A carga total é negativa: {carga_total:.4f}")
else:
    print("A carga total é neutra.")
