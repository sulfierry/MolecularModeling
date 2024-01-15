import matplotlib.pyplot as plt

def ler_cargas_mol2(arquivo_mol2):
    """ Lê as cargas parciais de um arquivo .mol2 """
    cargas = []
    with open(arquivo_mol2, 'r') as file:
        inicio = False
        for linha in file:
            if linha.startswith("@<TRIPOS>ATOM"):
                inicio = True
                continue
            elif linha.startswith("@<TRIPOS>"):
                inicio = False
            if inicio:
                partes = linha.split()
                if len(partes) > 8:
                    carga = float(partes[8])
                    cargas.append(carga)
    return cargas



def ajustar_cargas_hidrogenios(cargas, diferenca):
    """ Ajusta as cargas dos hidrogênios para corrigir a diferença com precisão de 6 casas decimais """
    hidrogenios = [i for i, carga in enumerate(cargas) if carga > 0]
    ajuste_por_hidrogenio = diferenca / len(hidrogenios)
    
    for i in hidrogenios:
        cargas[i] -= ajuste_por_hidrogenio
        # Arredondar cada carga para 6 casas decimais
        cargas[i] = round(cargas[i], 6)

    # Ajuste final para garantir que a soma seja 0.000000
    carga_total_ajustada = round(sum(cargas), 6)
    if carga_total_ajustada != 0.000000:
        cargas[hidrogenios[-1]] -= carga_total_ajustada

    return cargas
