import torch
import esm
import json
import numpy as np

class WeightExtract:
    def __init__(self, model_name: str):
        """
        Inicializa a classe carregando o modelo ESM correspondente.

        :param model_name: Nome do modelo ESM a ser carregado
        """
        self.model_name = model_name
        self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        self.weights_info = {}  # Dicionário para rastrear os pesos
        self.weights_data = {}  # Dicionário para armazenar os valores dos pesos

    def get_model_structure(self):
        """
        Exibe a estrutura da rede ESM, incluindo número de camadas e tipos de parâmetros.
        """
        print("\n🔎 Estrutura da Rede Neural:")
        print(self.model)

    def extract_weights(self):
        """
        Extrai os pesos da rede neural, incluindo seus valores exatos e dimensões.
        """
        print("\n📌 Extraindo pesos da rede...")

        for name, param in self.model.named_parameters():
            self.weights_info[name] = param.shape  # Armazena o formato dos pesos
            self.weights_data[name] = param.detach().cpu().numpy()  # Converte para NumPy para rastreamento
        
        print("✅ Extração de pesos concluída!")

    def display_weights(self):
        """
        Exibe a lista completa dos pesos extraídos, com nomes das camadas e suas dimensões.
        """
        print("\n📌 Mapeamento dos Pesos por Camada:")
        for layer, shape in self.weights_info.items():
            print(f"🔹 {layer}: {shape}")

    def save_weights(self, output_file: str):
        """
        Salva os pesos extraídos no formato `.pth` e `.json` para rastreabilidade.

        :param output_file: Nome base do arquivo onde os pesos serão salvos
        """
        # Salvar o estado bruto do modelo
        torch.save(self.model.state_dict(), f"{output_file}.pth")
        
        # Salvar os metadados dos pesos em JSON
        with open(f"{output_file}_metadata.json", "w") as f:
            json.dump({k: list(v) for k, v in self.weights_info.items()}, f, indent=4)

        # Salvar os pesos reais como NumPy para rastreamento exato
        np.savez(f"{output_file}_weights.npz", **self.weights_data)

        print(f"\n✅ Pesos salvos com rastreabilidade em:")
        print(f"  - '{output_file}.pth' (estado bruto da rede)")
        print(f"  - '{output_file}_metadata.json' (estrutura das camadas)")
        print(f"  - '{output_file}_weights.npz' (valores exatos dos pesos)")

    def trace_weight(self, layer_name: str, index: tuple):
        """
        Permite rastrear um peso específico dentro de uma camada.

        :param layer_name: Nome da camada (ex: 'encoder.layer.0.attention.self.query.weight')
        :param index: Índice do peso dentro do tensor (ex: (10, 20))
        """
        if layer_name not in self.weights_data:
            print(f"❌ Camada '{layer_name}' não encontrada!")
            return None
        
        try:
            value = self.weights_data[layer_name][index]
            print(f"\n🔍 Peso localizado: {layer_name}{index} = {value}")
            return value
        except IndexError:
            print(f"❌ Índice {index} fora do intervalo da matriz {self.weights_info[layer_name]}")
            return None

# =====================================================
#                   EXECUÇÃO NO MAIN
# =====================================================

def main():
    # 📌 Defina aqui o nome do modelo ESM desejado
    model_name = "esm2_t36_3B_UR50D"  # 🔧 Ajuste conforme necessário

    # 📌 Defina o nome do arquivo base para salvar os pesos
    output_file = "esm_model"

    # Criar instância da classe e executar métodos
    extractor = WeightExtract(model_name)

    extractor.get_model_structure()  # Exibir a estrutura do modelo
    extractor.extract_weights()  # Extrair os pesos
    extractor.display_weights()  # Exibir os pesos extraídos
    extractor.save_weights(output_file)  # Salvar os pesos

    # 📌 Exemplo de rastreabilidade: buscar um peso exato
    # layer = "encoder.layer.0.attention.self.query.weight"
    # index = (10, 20)  # Ajuste conforme necessário
    # extractor.trace_weight(layer, index)  # Buscar peso específico

if __name__ == "__main__":
    main()
