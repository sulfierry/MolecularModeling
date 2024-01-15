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


def formatar_linha_atom(linha, carga):
    """ Formata a linha do átomo com espaçamento adequado para a carga, mantendo 6 casas decimais """
    partes = linha.split()
    
    # Garante que -0.000000 seja formatado como 0.000000
    carga_formatada = f"{carga:.6f}" if carga != -0.000000 else "0.000000"
    partes[8] = carga_formatada
    return "{:<7s}{:<11s}{:>10s}{:>10s}{:>10s} {:<5s}{:>4s} {:<4s}{:>10s}\n".format(*partes)

def salvar_mol2_com_cargas_ajustadas(arquivo_original, arquivo_salvo, cargas_ajustadas):
    """ Salva um novo arquivo .mol2 com as cargas ajustadas """
    with open(arquivo_original, 'r') as file:
        linhas = file.readlines()

    indice_carga = 0
    for i in range(len(linhas)):
        if linhas[i].startswith("@<TRIPOS>ATOM"):
            inicio = i + 1
            break

    for i in range(inicio, len(linhas)):
        if linhas[i].startswith("@<TRIPOS>"):
            fim = i
            break
        partes = linhas[i].split()
        if len(partes) > 8:  # Se a linha tem uma carga
            linhas[i] = formatar_linha_atom(linhas[i], cargas_ajustadas[indice_carga])
            indice_carga += 1

    with open(arquivo_salvo, 'w') as file:
        file.writelines(linhas)


def armazenar_dados_cargas(arquivo_mol2):
    """ Armazena os dados dos átomos e suas cargas parciais de um arquivo .mol2 """
    dados_cargas = []
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
                    atomo = partes[1]
                    carga = float(partes[8])
                    dados_cargas.append((atomo, carga))
    return dados_cargas

def calcular_diferencas_estatisticas(dados_iniciais, dados_corrigidos):
    """ Calcula a diferença estatística para cada linha das colunas de cargas """
    diferenças = []
    for (atomo_ini, carga_ini), (atomo_cor, carga_cor) in zip(dados_iniciais, dados_corrigidos):
        diferença = carga_cor - carga_ini
        diferenças.append((atomo_ini, diferença))
    return diferenças


def main():
    input_name = 'ligand.mol2'
    output_name = 'ligand_adjusted.mol2'

    # Lendo e ajustando as cargas
    cargas_antes = ler_cargas_mol2(input_name)
    diferenca = sum(cargas_antes)
    cargas_ajustadas = ajustar_cargas_hidrogenios(cargas_antes.copy(), diferenca)

    # Calculando a diferença percentual
    soma_antes = sum(cargas_antes)
    soma_ajustadas = sum(cargas_ajustadas)
    diferenca_percentual = abs((soma_ajustadas - soma_antes) / soma_antes) * 100 if soma_antes != 0 else 0

    # Salvando o arquivo ajustado
    salvar_mol2_com_cargas_ajustadas(input_name, output_name, cargas_ajustadas)

    # Calculando e imprimindo as diferenças estatísticas
    dados_iniciais = armazenar_dados_cargas(input_name)
    dados_corrigidos = armazenar_dados_cargas(output_name)
    diferencas = calcular_diferencas_estatisticas(dados_iniciais, dados_corrigidos)

    for atomo, diferenca in diferencas:
        print(f"Atomo: {atomo}, Diferença: {diferenca:.6f}")

    # Plotando histogramas
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.hist(cargas_antes, bins=20, color='lightblue', alpha=0.6)
    plt.title("Antes do Ajuste")
    plt.xlabel("Carga Parcial")
    plt.ylabel("Frequência")

    plt.subplot(1, 2, 2)
    plt.hist(cargas_ajustadas, bins=20, color='salmon', alpha=0.6)
    plt.title("Após o Ajuste")
    plt.xlabel("Carga Parcial")
    plt.ylabel("Frequência")

    plt.suptitle("Histograma de Cargas Parciais")
    plt.figtext(0.5, 0.01, f'Diferença percentual: {diferenca_percentual:.2f}%', ha='center')
    plt.show()
    
    print(f"Diferença antes do ajuste: {soma_antes:.6f}")
    print(f"Diferença após o ajuste: {soma_ajustadas:.6f}")


# Execute a função main
if __name__ == "__main__":
    main()
