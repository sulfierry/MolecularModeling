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
