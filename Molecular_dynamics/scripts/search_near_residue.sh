#!/bin/bash

main() {
    # Define o caminho absoluto onde as pastas estão localizadas e outras variáveis
    export ROOT_DIR="/media/leon/FEDF-FDB3/way_analogues/01_gmmsb/run_15/water_remov"
    export base_name="gmmsb01_run15_wr" # os arquivos precisar ter _{i} apos o $base_name
    export contact_threshold_percentage="10"
    export ligand_distance="3.00"
    export ligand_name="LIG"
    export start_folder="1"
    export end_folder="10"

    # Chama a função que executa o protocolo
    execute_protocol
}

execute_protocol() {

    # Gera o script TCL antes de iniciar o loop
    generate_tcl_script
    
    echo "Iniciando protocolo..."

    # Loop pelas pastas
    for ((i=start_folder; i<=end_folder; i++)); do
        echo "Processando pasta $i..."

        # Verifica se a pasta existe antes de tentar acessá-la
        if [ -d "$ROOT_DIR/$i" ]; then
            # Acessa a pasta
            cd "$ROOT_DIR/$i"
            
            echo "Diretório atual: $(pwd)"
            # ls -l

            # Define os nomes dos arquivos de topologia e trajetória baseados na pasta atual
            TOPOLOGY_FILENAME="$base_name"_"${i}.prmtop"
            TRAJECTORY_FILENAME="$base_name"_"${i}.dcd"

            # Executa o VMD com o script TCL, passando os nomes dos arquivos como argumentos, se necessário
            vmd -dispdev text -e ../search_near_residue.tcl $TOPOLOGY_FILENAME $TRAJECTORY_FILENAME
            
            # Retorna ao diretório raiz antes de continuar o loop
            cd "$ROOT_DIR"
        else
            echo "A pasta $ROOT_DIR/$i não existe."
        fi
    done

    echo "Processamento completo."
    # Remove o script TCL após o uso
    rm search_near_residue.tcl
    
    near_residue_prevalence
    # Executa o script Python para análise de prevalência dos resíduos
    if [ -f "near_residue_prevalence.py" ]; then
        python near_residue_prevalence.py
        # Remove o script Python após o uso, se necessário
        rm near_residue_prevalence.py
    else
        echo "O arquivo near_residue_prevalence.py não foi encontrado."
    fi
}




# Função para gerar o script search_near_residue.tcl
generate_tcl_script() {
    cat << EOF > search_near_residue.tcl
# vmd -dispdev text -e script.tcl 
# Find the number of atoms within 2-3 A of the ligand. 

# Define o nome do ligante
set ligand_name $ligand_name

# Cria uma seleção de todos os resíduos de proteína próximos ao ligante
set closer_residue [atomselect 0 "protein within $ligand_distance of resname \$ligand_name"]

# Abre o arquivo para escrita
set outfile [open "residues_near_ligand.dat" w]

# Obtém o número total de quadros na trajetória
set num_frames [molinfo 0 get numframes]

# Itera sobre cada quadro da simulação
for {set frame_number 0} {\$frame_number < \$num_frames} {incr frame_number} {
    # Atualiza a seleção para o quadro atual
    \$closer_residue frame \$frame_number
    \$closer_residue update

    # Obtém os resíduos únicos próximos ao ligante neste quadro
    set residues [\$closer_residue get {resname resid}]
    set unique_residues [lsort -unique \$residues]

    # Escreve o número do quadro e o número de contatos únicos
    puts \$outfile "Frame \$frame_number: [llength \$unique_residues] unique residues near \$ligand_name"

    # Para cada resíduo único, imprime sua identificação
    foreach {resname resid} \$unique_residues {
        puts \$outfile "    Residue \$resname \$resid"
    }
}

# Fecha o arquivo
close \$outfile

# Encerra o script VMD
quit
EOF
}


near_residue_prevalence() {
    cat << EOF > near_residue_prevalence.py
import os
import re
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict

class ResiduesPrevalence:

    def __init__(self, ligand_name, threshold=10, color='blue', root_dir='.'):
        self.ligand_name = ligand_name
        self.file_path = os.path.join(root_dir, 'residues_near_ligand_total.dat')
        self.threshold = threshold
        self.color = color
        self.root_dir = root_dir

    def merge_data_files(self):
        """Realiza o merge dos arquivos residues_near_ligand.dat em todas as pastas."""
        unique_residues_by_frame = defaultdict(set)

        # Substitua [1, 11] pelo intervalo correto de pastas, se necessário
        for i in range(${start_folder}, ${end_folder}+1):
            folder_path = os.path.join(self.root_dir, str(i))
            file_path = os.path.join(folder_path, "residues_near_ligand.dat")

            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    current_frame = None
                    for line in file:
                        if line.startswith("Frame"):
                            current_frame = line.split()[1].strip(':')
                        else:
                            residue_info = line.strip().split()[1:] # Remove 'Residue' e apenas extrai os resíduos
                            for j in range(0, len(residue_info), 2):
                                residue = f"{residue_info[j]} {residue_info[j+1]}"
                                unique_residues_by_frame[current_frame].add(residue)

        # Salvar os dados consolidados em um novo arquivo
        with open(self.file_path, 'w') as output_file:
            for frame, residues in sorted(unique_residues_by_frame.items(), key=lambda x: int(x[0])):
                output_file.write(f"Frame {frame}: {len(residues)} unique residues near ligand\n")
                for residue in sorted(residues):
                    output_file.write(f"    Residue {residue}\n")


    def process_data(self):
        """Processa o arquivo de dados para calcular a prevalência de cada resíduo."""
        residue_count = defaultdict(int)
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('    Residue '):
                    residues = re.findall(r'([A-Z]+ \d+)', line)
                    for residue in residues:
                        residue_count[residue] += 1

        total_frames = max(int(line.split()[1].rstrip(':')) for line in open(self.file_path) if line.startswith("Frame")) + 1
        self.residue_prevalence = {residue: count / total_frames for residue, count in residue_count.items()}

    def filter_by_threshold(self):
        """Filtra os resíduos com prevalência acima do threshold especificado."""
        self.filtered_data = {res: prevalence for res, prevalence in self.residue_prevalence.items() if prevalence * 100 >= self.threshold}

    def plot_data(self):
        """Gera o gráfico com os dados processados e filtrados, ordenados pelo número do resíduo."""
        sorted_data = sorted(self.filtered_data.items(), key=lambda x: int(x[0].split()[1]))
        labels = [f"{res[0].split()[0]} {res[0].split()[1]}" for res in sorted_data]
        values = [prevalence * 100 for _, prevalence in sorted_data]
        plt.figure(figsize=(10, 8))
        plt.bar(labels, values, color=self.color)
        plt.xlabel('Amino Acid and Number')
        plt.ylabel('Prevalence (%)')
        plt.xticks(rotation=45, ha="right")
        plt.title('Amino Acid Prevalence Near Ligand')
        plt.tight_layout()
        plt.savefig('aac_prevalence.png')
        plt.show()


    @classmethod
    def main(cls):
        ligand_name = "${ligand_name}"
        threshold = ${contact_threshold_percentage}
        color = 'royalblue'
        root_dir = '.'  # Caminho para o diretório raiz onde estão localizadas as pastas

        # Cria uma instância da classe com as configurações desejadas
        instance = cls(ligand_name, threshold, color, root_dir)

        # Realiza o merge dos arquivos de dados
        instance.merge_data_files()

        # Processa o arquivo de dados gerado
        instance.process_data()

        # Filtra os dados com base no threshold especificado
        instance.filter_by_threshold()

        # Plota os dados filtrados
        instance.plot_data()

if __name__ == '__main__':
    ResiduesPrevalence.main()

EOF
}

main
