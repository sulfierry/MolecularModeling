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
