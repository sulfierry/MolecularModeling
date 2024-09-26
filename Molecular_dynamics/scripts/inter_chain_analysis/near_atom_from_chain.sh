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
    python_script=$(generate_python_script $TRESHOLD_PREVALENCE_INTERACTION)

    # Salva o script Python gerado em um arquivo temporário
    python_script_file="${ROOT_DIR}/process_interaction.py"
    echo "$python_script" > "$python_script_file"

    # Executa o script Python gerado
    python3 "$python_script_file"

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
import glob
import pandas as pd
import matplotlib.pyplot as plt

class ProcessInteractions:

    def __init__(self, input_dir, output_file):
        self.input_dir = input_dir
        self.output_file = output_file

    def reformat_data(self):
        dat_files = sorted(glob.glob(os.path.join(self.input_dir, 'detailed_interactions_frame_*.dat')))
        with open(self.output_file, 'w') as outfile:
            outfile.write("frame\tchainA_residue\tchainA_residue_number\tchainA_atom\tchainB_residue\tchainB_residue_number\tchainB_atom\tdistance\n")
            for dat_file in dat_files:
                with open(dat_file, 'r') as infile:
                    for line in infile:
                        parts = line.strip().split()
                        frame = parts[1].strip(':')
                        chainA_residue = parts[3]
                        chainA_residue_number = parts[4]
                        chainA_atom = parts[5]
                        chainB_residue = parts[8]
                        chainB_residue_number = parts[9]
                        chainB_atom = parts[10]
                        distance = line.split('Distance: ')[1].strip()
                        output_line = f"{frame}\t{chainA_residue}\t{chainA_residue_number}\t{chainA_atom}\t{chainB_residue}\t{chainB_residue_number}\t{chainB_atom}\t{distance}\n"
                        outfile.write(output_line)
        print(f"Processamento concluído. Todos os dados foram salvos em: {self.output_file}\n")

    def calculate_prevalence(self, prevalence_threshold=50):
        data = pd.read_csv(self.output_file, sep='\t')
        interaction_counts = data.groupby(['chainA_residue', 'chainA_residue_number', 'chainA_atom', 'chainB_residue', 'chainB_residue_number', 'chainB_atom']).size()
        total_frames = data['frame'].nunique()
        prevalence = (interaction_counts / total_frames) * 100
        prevalence_filtered = prevalence[prevalence >= prevalence_threshold].reset_index(name='prevalence')
        self.prevalence_data = prevalence_filtered
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print("-------------------------------------------------------------------------------------------------------------------")
            print(self.prevalence_data)
            print("-------------------------------------------------------------------------------------------------------------------\n")

    def plot_prevalence(self):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return
        self.prevalence_data['chainB_residue_number'] = pd.to_numeric(self.prevalence_data['chainB_residue_number'])
        self.prevalence_data.sort_values(by='chainB_residue_number', inplace=True)
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
        plt.show()

    def identify_rotatable_atoms(self, prevalence_threshold=50):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return
        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= prevalence_threshold]
        print(f"Rotatable atoms: \n{rotatable_atoms}")
        print("-------------------------------------------------------------------------------------------------------------------\n")

    @staticmethod
    def main():
        input_dir = 'vmd_frames'
        output_file = 'detailed_interaction.tsv'
        prevalence_threshold = 80
        process_interactions = ProcessInteractions(input_dir, output_file)
        process_interactions.reformat_data()
        process_interactions.calculate_prevalence(prevalence_threshold)
        process_interactions.plot_prevalence()
        process_interactions.identify_rotatable_atoms(prevalence_threshold)

if __name__ == '__main__':
    ProcessInteractions.main()
EOF
}

# Verifica se o script é o principal
if [ "$0" = "${BASH_SOURCE[0]}" ]; then
    main "$@"
fi
