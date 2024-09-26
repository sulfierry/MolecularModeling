#!/bin/bash

# Este script deve ser ajustado de acodo com as novas classes para processamente e visualziacao das interacoes


# Função principal
main() {
    # Configuração inicial
    ROOT_DIR="$(pwd)"
    BASE_NAME="top1_wr_1"
    TOPOLOGY="$BASE_NAME.prmtop"
    COORDINATES="0_1_output.dcd"
    TOTAL_FRAMES="3224"
    SLICE="10"
    CHAIN_A_RESIDUES="resid 1 to 1323"
    CHAIN_B_RESIDUES="resid 1324 to 1905"
    CUTOFF_DISTANCE="3.0"
    TRESHOLD_PREVALENCE_INTERACTION="80"

    process_frames
}

process_frames() {
    # Processa os frames em blocos
    for ((start=0; start<TOTAL_FRAMES; start+=SLICE)); do
        end=$((start + SLICE - 1))
        if [ $end -ge $TOTAL_FRAMES ]; then
            end=$((TOTAL_FRAMES - 1))
        fi

        # Gera o script TCL para o bloco atual de frames
        tcl_script=$(generate_tcl_script $start $end)

        # Salva o script TCL gerado em um arquivo temporário
        tcl_script_file="${ROOT_DIR}/temp_tcl_script.tcl"
        echo "$tcl_script" > "$tcl_script_file"

        # Executa o VMD com o script TCL gerado
        vmd -dispdev text -e "$tcl_script_file" "$TOPOLOGY" "$COORDINATES"

        # Limpa o arquivo temporário do script TCL
        rm "$tcl_script_file"
    done

    # Gera o script Python dinamicamente
    python_script=$(generate_python_script)

    # Salva o script Python gerado em um arquivo temporário
    python_script_file="${ROOT_DIR}/process_interaction.py"
    echo "$python_script" > "$python_script_file"

    # Executa o script Python gerado
    python3 "$python_script_file" --input_dir vmd_frames --output_file detailed_interactions.tsv --prevalence_threshold $TRESHOLD_PREVALENCE_INTERACTION --top 30

    # Limpa o arquivo temporário do script Python
    rm "$python_script_file"
}

# Função para gerar o script TCL dinamicamente
generate_tcl_script() {
    local start_frame=$1
    local end_frame=$2
    cat << EOF
# Verifica e cria a pasta "vmd_frames" se ela não existir
if {![file exists vmd_frames]} {
    file mkdir vmd_frames
}

# Defina a distância de corte
set cutoff_distance ${CUTOFF_DISTANCE}

# Selecione os resíduos das cadeias A e B
set chainA_residues "${CHAIN_A_RESIDUES}"
set chainB_residues "${CHAIN_B_RESIDUES}"

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação dentro do bloco definido
for {set frame_number $start_frame} {\$frame_number <= $end_frame} {incr frame_number} {
    # Atualiza para o quadro atual
    molinfo 0 set frame \$frame_number

    # Abre o arquivo para escrita para o quadro atual dentro da pasta "vmd_frames"
    set outfile [open "vmd_frames/detailed_interactions_frame_\$frame_number.dat" w]

    # Cria seleções para as cadeias A e B
    set chainA_atoms [atomselect 0 "\$chainA_residues"]
    set chainB_atoms [atomselect 0 "\$chainB_residues"]

    # Itera sobre os átomos da cadeia A
    foreach chainA_atom [\$chainA_atoms get {index name resname resid}] {
        set chainA_index [lindex \$chainA_atom 0]
        set chainA_atomname [lindex \$chainA_atom 1]
        set chainA_resname [lindex \$chainA_atom 2]
        set chainA_resid [lindex \$chainA_atom 3]
        
        # Cria uma seleção de átomos da cadeia B próximos ao átomo atual da cadeia A
        set nearby_atoms [atomselect 0 "same residue as (within \$cutoff_distance of index \$chainA_index) and \$chainB_residues"]
        set nearby_atoms_info [\$nearby_atoms get {index element resname resid name}]
        
        foreach nearby_atom \$nearby_atoms_info {
            set nearby_index [lindex \$nearby_atom 0]
            set nearby_element [lindex \$nearby_atom 1]
            set nearby_resname [lindex \$nearby_atom 2]
            set nearby_resid [lindex \$nearby_atom 3]
            set nearby_atomname [lindex \$nearby_atom 4]
            
            # Calcula a distância entre o átomo da cadeia A e o átomo da cadeia B
            set distance [measure bond [list \$chainA_index \$nearby_index]]

            # Escreve as informações no arquivo, incluindo a distância
            puts \$outfile "Frame \$frame_number: Chain A \$chainA_resname \$chainA_resid \$chainA_atomname interacts with Chain B \$nearby_resname \$nearby_resid \$nearby_atomname Distance: \$distance"
        }
    }

    # Fecha o arquivo
    close \$outfile
}

# Encerra o script VMD
quit
EOF
}

generate_python_script() {
    cat << EOF
import os
import sys
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

class ProcessInteractions:
    def __init__(self, input_dir, output_file, prevalence_threshold):
        self.input_dir = input_dir
        self.output_file = output_file
        self.prevalence_threshold = prevalence_threshold
        self.data = None
        self.prevalence_data = None
        self.filtered_data = None

    def load_data(self):
        dat_files = sorted(glob.glob(os.path.join(self.input_dir, 'detailed_interactions_frame_*.dat')))
        print(f"Número de arquivos .dat encontrados: {len(dat_files)}")
        
        data_list = []
        for dat_file in dat_files:
            with open(dat_file, 'r') as infile:
                for line in infile:
                    parts = line.strip().split()
                    frame = int(parts[1].strip(':'))
                    data_list.append({
                        'frame': frame,
                        'chainA_residue': parts[4],
                        'chainA_residue_number': int(parts[5]),
                        'chainA_atom': parts[6],
                        'chainB_residue': parts[11],
                        'chainB_residue_number': int(parts[12]),
                        'chainB_atom': parts[13],
                        'distance': float(parts[-1])
                    })
        
        self.data = pd.DataFrame(data_list)
        print(f"Total de interações carregadas: {len(self.data)}")
        print("Átomos únicos da cadeia B nos dados carregados:")
        print(self.data['chainB_atom'].unique())

    def calculate_prevalence(self):
        if self.data is None:
            print("Dados não carregados. Execute load_data() primeiro.")
            return

        print(f"Dimensões do DataFrame: {self.data.shape}")
        print("Primeiras linhas dos dados:")
        print(self.data.head())
        
        total_frames = self.data['frame'].nunique()
        print(f"\nNúmero total de frames: {total_frames}")
        
        interaction_counts = self.data.groupby(['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom']).size()
        prevalence = (interaction_counts / total_frames) * 100
        
        self.prevalence_data = prevalence[prevalence >= self.prevalence_threshold].reset_index(name='prevalence')
        
        print(f"\nNúmero de interações que atingiram o limiar de prevalência ({self.prevalence_threshold}%): {len(self.prevalence_data)}")
        
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print("\nInterações que atingiram o limiar de prevalência:")
            print(self.prevalence_data)

        print("\nÁtomos únicos da cadeia B nos dados de prevalência:")
        print(self.prevalence_data['chainB_atom'].unique())

        # Filtrar os dados originais com base nas interações que atingiram o limiar de prevalência
        self.filtered_data = self.data.merge(self.prevalence_data, on=['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom'])
        
        print("\nÁtomos únicos da cadeia B nos dados filtrados:")
        print(self.filtered_data['chainB_atom'].unique())

        # Verificação explícita dos átomos da cadeia B
        unique_chainB_atoms = self.filtered_data['chainB_atom'].unique()
        print(f"Átomos únicos da cadeia B nos dados filtrados: {unique_chainB_atoms}")
        if len(unique_chainB_atoms) == 1 and unique_chainB_atoms[0] == 'N':
            print("AVISO: Apenas o átomo 'N' está presente para a cadeia B nos dados filtrados.")

        # Salvar os dados filtrados
        self.filtered_data.to_csv(self.output_file, sep='\t', index=False)
        print(f"\nDados filtrados salvos em: {self.output_file}")
        print(f"Número de linhas no arquivo de saída: {len(self.filtered_data)}")

    def plot_prevalence(self):
        if self.prevalence_data is None:
            print("Dados de prevalência não calculados. Por favor, execute calculate_prevalence primeiro.")
            return
        
        self.prevalence_data['chainB_residue_number'] = pd.to_numeric(self.prevalence_data['chainB_residue_number'], errors='coerce')
        self.prevalence_data = self.prevalence_data.sort_values(by='chainB_residue_number', na_position='last')
        
        self.prevalence_data['description'] = self.prevalence_data.apply(
            lambda x: f"{x['chainB_residue_number']}{x['chainB_residue']}({x['chainB_atom']}) - {x['chainA_residue_number']}{x['chainA_residue']}({x['chainA_atom']})", axis=1)
        
        plt.figure(figsize=(14, 10))
        plt.bar(self.prevalence_data['description'], self.prevalence_data['prevalence'], color='skyblue')
        plt.ylabel('Prevalência (%)')
        plt.xlabel('Contato Átomo da Cadeia B - Átomo da Cadeia A')
        plt.xticks(rotation=90, ha="right")
        plt.title('Prevalência de Contato Átomo da Cadeia A - Átomo da Cadeia B')
        plt.tight_layout()
        plt.savefig('prevalence.png')
        print("Gráfico de prevalência salvo como 'prevalence.png'")
        plt.show()

    def identify_rotatable_atoms(self):
        if self.prevalence_data is None:
            print("Dados de prevalência não calculados. Por favor, execute calculate_prevalence primeiro.")
            return
        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= self.prevalence_threshold]
        print(f"Átomos rotacionáveis (prevalência >= {self.prevalence_threshold}%):")
        print(rotatable_atoms)
        print(f"Total de átomos rotacionáveis identificados: {len(rotatable_atoms)}")

    def run(self):
        self.load_data()
        self.calculate_prevalence()
        self.identify_rotatable_atoms()
        # Uncomment the following line if you want to plot the prevalence
        # self.plot_prevalence()

class ViewInteractions:
    def __init__(self, input_file, top_interactions):
        self.input_file = input_file
        self.top_interactions = top_interactions
        self.data = None
        self.interaction_grouped = None

    def load_data(self):
        if not os.path.exists(self.input_file):
            print(f"Erro: O arquivo {self.input_file} não existe.")
            exit(1)
        try:
            self.data = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Erro ao ler o arquivo: {e}")
            exit(1)

    def process_data(self):
        print("Primeiras linhas do arquivo de entrada:")
        print(self.data.head())

        print("\nÁtomos únicos na cadeia A (antes do agrupamento):")
        print(self.data['chainA_atom'].unique())
        print("\nÁtomos únicos na cadeia B (antes do agrupamento):")
        print(self.data['chainB_atom'].unique())

        print("\nContagem de átomos na cadeia B:")
        print(self.data['chainB_atom'].value_counts())

        self.data['interaction'] = (self.data['chainA_residue'] + ' (' + self.data['chainA_residue_number'].astype(str) + ') - ' +
                                    self.data['chainB_residue'] + ' (' + self.data['chainB_residue_number'].astype(str) + ')')

        self.interaction_grouped = self.data.groupby('interaction').apply(self.aggregate_interaction_data).reset_index()
        self.interaction_grouped = self.interaction_grouped.sort_values(by='prevalence', ascending=False)

        print("\nPrimeiras interações mais prevalentes:")
        print(self.interaction_grouped.head())

        print("\nÁtomos únicos na cadeia A (após o agrupamento):")
        print(self.interaction_grouped['chainA_atom'].unique())
        print("\nÁtomos únicos na cadeia B (após o agrupamento):")
        print(self.interaction_grouped['chainB_atom'].unique())

        print("\nAmostra de linhas para chainB_atom:")
        print(self.interaction_grouped[['interaction', 'chainB_atom']].sample(min(10, len(self.interaction_grouped))))

    @staticmethod
    def aggregate_interaction_data(group):
        chainB_atoms = group['chainB_atom'].value_counts()
        most_common_chainB_atom = chainB_atoms.index[0] if not chainB_atoms.empty else 'Unknown'
        return pd.Series({
            'prevalence': group['prevalence'].mean(),
            'distance': group['distance'].mean(),
            'chainA_atom': group['chainA_atom'].mode().iloc[0],
            'chainB_atom': most_common_chainB_atom
        })

    def save_csv(self):
        output_csv = 'aminoacid_interactions_prevalence.csv'
        self.interaction_grouped.to_csv(output_csv, index=False)
        print(f"\nArquivo CSV salvo como '{output_csv}'")

    def plot_interactions(self):
        plt.figure(figsize=(12, 10))
        data_to_plot = self.interaction_grouped.head(self.top_interactions)
        sns.barplot(x='prevalence', y='interaction', data=data_to_plot, palette='viridis')
        plt.title(f'Top {len(data_to_plot)} Most Prevalent Interactions Between Chains')
        plt.xlabel('Average Prevalence (%)')
        plt.ylabel('Interaction (Residue-Residue)')
        plt.tight_layout()

        output_png = f'prev_interaction_top_{len(data_to_plot)}.png'
        plt.savefig(output_png, bbox_inches='tight')
        print(f"\nGráfico salvo como '{output_png}'")
        plt.show()

    def run(self):
        self.load_data()
        self.process_data()
        self.save_csv()
        self.plot_interactions()

def main():
    parser = argparse.ArgumentParser(description='Process and visualize protein chain interactions.')
    parser.add_argument('--input_dir', type=str, default='vmd_frames', help='Directory containing input .dat files')
    parser.add_argument('--output_file', type=str, default='detailed_interactions.tsv', help='Output file name')
    parser.add_argument('--prevalence_threshold', type=float, default=50, help='Prevalence threshold for interactions')
    parser.add_argument('--top', type=int, default=30, choices=[10, 20, 30, 40, 50], help='Number of top interactions to visualize')
    args = parser.parse_args()

    # Process interactions
    process_interactions = ProcessInteractions(args.input_dir, args.output_file, args.prevalence_threshold)
    process_interactions.run()

    # View interactions
    view_interactions = ViewInteractions(args.output_file, args.top)
    view_interactions.run()

if __name__ == '__main__':
    main()
EOF
}

# Verifica se o script é o principal
if [ "$0" = "${BASH_SOURCE[0]}" ]; then
    main "$@"
fi
