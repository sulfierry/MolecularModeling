import matplotlib.pyplot as plt

class ChargeAdjust:

    def __init__(self, input_name):
        self.input_name = input_name
        self.output_name = input_name.replace('.mol2', '_adjusted.mol2')
    
    def ler_cargas_mol2(self):
        cargas = []
        with open(self.input_name, 'r') as file:
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
    def formatar_linha_atom(self, linha, carga):
        partes = linha.split()
        carga_formatada = f"{carga:.6f}" if carga != -0.000000 else "0.000000"
        partes[8] = carga_formatada
        return "{:<7s}{:<11s}{:>10s}{:>10s}{:>10s} {:<5s}{:>4s} {:<4s}{:>10s}\n".format(*partes)

    def salvar_mol2_com_cargas_ajustadas(self, cargas_ajustadas):
        with open(self.input_name, 'r') as file:
            linhas = file.readlines()

        indice_carga = 0
        for i in range(len(linhas)):
            if linhas[i].startswith("@<TRIPOS>ATOM"):
                inicio = i + 1
                break

        for i in range(inicio, len(linhas)):
            if linhas[i].startswith("@<TRIPOS>"):
                break
            partes = linhas[i].split()
            if len(partes) > 8:
                linhas[i] = self.formatar_linha_atom(linhas[i], cargas_ajustadas[indice_carga])
                indice_carga += 1

        with open(self.output_name, 'w') as file:
            file.writelines(linhas)

    def armazenar_dados_cargas(self, arquivo_mol2):
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


    def calcular_diferencas_estatisticas(self, dados_iniciais, dados_corrigidos):
        diferenças = []
        for (atomo_ini, carga_ini), (atomo_cor, carga_cor) in zip(dados_iniciais, dados_corrigidos):
            diferença = carga_cor - carga_ini
            diferenças.append((atomo_ini, diferença))
        return diferenças

    def plotar_histogramas(self, cargas_antes, cargas_ajustadas):
        plt.figure(figsize=(12, 6))

        # Histograma antes do ajuste
        plt.subplot(1, 2, 1)
        plt.hist(cargas_antes, bins=20, color='lightblue', alpha=0.7)
        #plt.axvline(x=sum(cargas_antes)/len(cargas_antes), color='darkblue', linestyle='dashed', linewidth=1)
        plt.title("Antes do Ajuste")
        plt.xlabel("Carga Parcial")
        plt.ylabel("Frequência")

        # Histograma após o ajuste
        plt.subplot(1, 2, 2)
        plt.hist(cargas_ajustadas, bins=20, color='salmon', alpha=0.7)
        #plt.axvline(x=sum(cargas_ajustadas)/len(cargas_ajustadas), color='darkred', linestyle='dashed', linewidth=1)
        plt.title("Após o Ajuste")
        plt.xlabel("Carga Parcial")
        plt.ylabel("Frequência")

        plt.suptitle("Histograma de Cargas Parciais: Antes e Depois do Ajuste")
        plt.show()


    def main(self):
        cargas_antes = self.ler_cargas_mol2()
        soma_cargas_antes = sum(cargas_antes)
        print(f"Carga total antes do ajuste: {soma_cargas_antes:.6f}")

        cargas_ajustadas = self.ajustar_cargas_hidrogenios(cargas_antes)
        soma_cargas_ajustadas = sum(cargas_ajustadas)
        print(f"Carga total após o ajuste: {soma_cargas_ajustadas:.6f}")

        self.salvar_mol2_com_cargas_ajustadas(cargas_ajustadas)

        dados_iniciais = self.armazenar_dados_cargas(self.input_name)
        dados_corrigidos = self.armazenar_dados_cargas(self.output_name)
        diferencas = self.calcular_diferencas_estatisticas(dados_iniciais, dados_corrigidos)

        for atomo, diferenca in diferencas:
            print(f"Atomo: {atomo}, Diferença: {diferenca:.6f}")

        self.plotar_histogramas(cargas_antes, cargas_ajustadas)

if __name__ == "__main__":
    charge_adjuster = ChargeAdjust('ligand.mol2')
    charge_adjuster.main()
