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

    def ajustar_cargas_hidrogenios(self, cargas):
        diferenca = sum(cargas)
        hidrogenios = [i for i, carga in enumerate(cargas) if carga > 0]
        ajuste_total = round(diferenca) - diferenca
        ajuste_por_hidrogenio = ajuste_total / len(hidrogenios)
        
        for i in hidrogenios:
            cargas[i] += ajuste_por_hidrogenio
            cargas[i] = round(cargas[i], 6)

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


    def calcular_diferencas_percentuais(self, dados_iniciais, dados_corrigidos):
        diferenças_percentuais = []
        for (atomo_ini, carga_ini), (atomo_cor, carga_cor) in zip(dados_iniciais, dados_corrigidos):
            if 'H' in atomo_ini:  # Considerando apenas hidrogênios
                diferença_percentual = ((carga_cor - carga_ini) / carga_ini) * 100 if carga_ini != 0 else 0
                diferenças_percentuais.append((atomo_ini, diferença_percentual))
        return diferenças_percentuais

    def plotar_diferencas_percentuais(self, diferencas_percentuais):
        atomos, diferencas = zip(*diferencas_percentuais)
        plt.figure(figsize=(12, 6))
        plt.bar(atomos, diferencas, color='purple')
        plt.xlabel("Átomos de Hidrogênio")
        plt.ylabel("Diferença Percentual (%)")
        plt.title("Diferença Percentual nas Cargas Parciais dos Hidrogênios")
        plt.xticks(rotation=45)
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
        diferencas_percentuais = self.calcular_diferencas_percentuais(dados_iniciais, dados_corrigidos)
        self.plotar_diferencas_percentuais(diferencas_percentuais)


if __name__ == "__main__":
    charge_adjuster = ChargeAdjust('ligand.mol2')
    charge_adjuster.main()
