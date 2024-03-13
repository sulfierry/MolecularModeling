#!/bin/bash

# Função principal
main() {
    # Configuração inicial
    ROOT_DIR="$(pwd)"
    BASE_NAME="gmmsb01_run15_wr_10"
    TOPOLOGY="$BASE_NAME.prmtop"
    COORDINATES="$BASE_NAME.dcd"
    TOTAL_FRAMES=1000
    SLICE=500
    LIGAND_NAME="LIG"
    CUTOFF_DISTANCE="3.0"
    TRESHOLD_PREVALENCE_INTERACTION="50"

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

# Defina o nome do ligante
set ligand_name "${LIGAND_NAME}"

# Defina a distância de corte
set cutoff_distance ${CUTOFF_DISTANCE}

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação dentro do bloco definido
for {set frame_number $start_frame} {\$frame_number <= $end_frame} {incr frame_number} {
    # Atualiza para o quadro atual
    molinfo 0 set frame \$frame_number

    # Abre o arquivo para escrita para o quadro atual dentro da pasta "vmd_frames"
    set outfile [open "vmd_frames/detailed_interactions_frame_\$frame_number.dat" w]

    # Cria uma seleção de todos os átomos do ligante
    set ligand_atoms [atomselect 0 "resname \$ligand_name"]

    # Obtém informações detalhadas dos átomos do ligante
    foreach ligand_atom [\$ligand_atoms get {index name}] {
        set ligand_index [lindex \$ligand_atom 0]
        set ligand_atomname [lindex \$ligand_atom 1]
        
        # Cria uma seleção de todos os átomos de cadeias laterais de aminoácidos próximos a este átomo do ligante
        set sidechain_atoms [atomselect 0 "sidechain within \$cutoff_distance of index \$ligand_index"]
        set sidechain_atoms_info [\$sidechain_atoms get {index element resname resid name}]
        
        foreach sidechain_atom \$sidechain_atoms_info {
            set sidechain_index [lindex \$sidechain_atom 0]
            set sidechain_element [lindex \$sidechain_atom 1]
            set sidechain_resname [lindex \$sidechain_atom 2]
            set sidechain_resid [lindex \$sidechain_atom 3]
            set sidechain_atomname [lindex \$sidechain_atom 4]
            
            # Calcula a distância entre o átomo do ligante e o átomo do aminoácido
            set distance [measure bond [list \$ligand_index \$sidechain_index]]

            # Escreve as informações no arquivo, incluindo a distância
            puts \$outfile "Frame \$frame_number: Ligand atom \$ligand_atomname index \$ligand_index interacts with \$sidechain_resname \$sidechain_resid \$sidechain_atomname Distance: \$distance"
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
            outfile.write("frame\tligand_atom\tligand_atom_index\taminoacid\taminoacid_atom\taminoacid_number\tdistance\n")
            for dat_file in dat_files:
                with open(dat_file, 'r') as infile:
                    for line in infile:
                        parts = line.strip().split()
                        frame = parts[1].strip(':')
                        ligand_atom = parts[4]
                        ligand_atom_index = parts[6]
                        aminoacid = parts[9]
                        aminoacid_number = parts[10]
                        aminoacid_atom = parts[11]
                        distance = line.split('Distance: ')[1].strip()
                        output_line = f"{frame}\t{ligand_atom}\t{ligand_atom_index}\t{aminoacid}\t{aminoacid_atom}\t{aminoacid_number}\t{distance}\n"
                        outfile.write(output_line)
        print(f"Processamento concluído. Todos os dados foram salvos em: {self.output_file}")

    def calculate_prevalence(self, prevalence_threshold=50):
        data = pd.read_csv(self.output_file, sep='\t')
        interaction_counts = data.groupby(['ligand_atom', 'ligand_atom_index', 'aminoacid', 'aminoacid_atom', 'aminoacid_number']).size()
        total_frames = data['frame'].nunique()
        prevalence = (interaction_counts / total_frames) * 100
        prevalence_filtered = prevalence[prevalence >= prevalence_threshold].reset_index(name='prevalence')
        self.prevalence_data = prevalence_filtered
        if self.prevalence_data.empty:
            print("Nenhuma interação atingiu o limiar de prevalência especificado.")
        else:
            print(self.prevalence_data)

    def plot_prevalence(self):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return
        self.prevalence_data['aminoacid_number'] = pd.to_numeric(self.prevalence_data['aminoacid_number'])
        self.prevalence_data.sort_values(by='aminoacid_number', inplace=True)
        self.prevalence_data['description'] = self.prevalence_data.apply(
            lambda x: f"{x['aminoacid_number']}{x['aminoacid']}({x['aminoacid_atom']}) - {x['ligand_atom']}(LIG)", axis=1)
        plt.figure(figsize=(14, 10))
        plt.bar(self.prevalence_data['description'], self.prevalence_data['prevalence'], color='skyblue')
        plt.ylabel('Prevalência (%)')
        plt.xlabel('Contato Átomo do Aminoácido - Átomo do Ligante')
        plt.xticks(rotation=90, ha="right")
        plt.title('Prevalência de Contato Átomo do Ligante - Átomo do Aminoácido')
        plt.tight_layout()
        plt.savefig('prevalence.png')
        plt.show()

    def identify_rotatable_atoms(self, prevalence_threshold=50):
        if not hasattr(self, 'prevalence_data'):
            print("Prevalence data not calculated. Please run calculate_prevalence first.")
            return
        rotatable_atoms = self.prevalence_data[self.prevalence_data['prevalence'] >= prevalence_threshold]
        print(f"Rotatable atoms: \n{rotatable_atoms}")

    @staticmethod
    def main():
        input_dir = 'vmd_frames'
        output_file = 'detailed_interaction.tsv'
        prevalence_threshold = $TRESHOLD_PREVALENCE_INTERACTION
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
