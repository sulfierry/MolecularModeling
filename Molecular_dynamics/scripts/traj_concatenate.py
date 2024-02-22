import MDAnalysis as mda
from MDAnalysis.analysis import align
import os
import shutil
from MDAnalysis.coordinates.DCD import DCDWriter

def concatenate_and_align_trajectories(base_path, output_path):
    # Garantir a criação do diretório de saída
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Caminho para a topologia inicial
    topology_path = os.path.join(base_path, "1/5cc8_wr_1.prmtop")
    # Lista de trajetórias a serem carregadas e concatenadas
    trajectory_paths = [os.path.join(base_path, f"{i}/5cc8_wr_{i}.dcd") for i in range(1, 11)]
    
    # Carregar a primeira trajetória para inicializar o universo
    u = mda.Universe(topology_path, trajectory_paths[0])
    
    # Selecionar átomos para alinhamento, ajuste a seleção conforme necessário
    ref_atoms = u.select_atoms('backbone')  # Exemplo: seleção da espinha dorsal

    # Criar um universo de referência com a estrutura do frame inicial
    ref_universe = mda.Universe(topology_path, trajectory_paths[0])
    ref_frame = ref_universe.trajectory[0]
    ref_positions = ref_universe.select_atoms('backbone').positions

    # Preparar o writer para salvar a trajetória alinhada
    with DCDWriter(os.path.join(output_path, "all_traj_aligned.dcd"), n_atoms=u.atoms.n_atoms) as W:
        for traj_path in trajectory_paths:
            # Atualizar o universo para a trajetória atual
            u.load_new(traj_path)
            for ts in u.trajectory:
                # Alinhar ao frame de referência
                align.alignto(u, ref_universe, select='backbone')
                # Escrever o frame alinhado
                W.write(u)
    
    # Copiar a topologia para a pasta de saída
    shutil.copy(topology_path, os.path.join(output_path, "5cc8_wr_1.prmtop"))

if __name__ == "__main__":
    base_path = "./"  # Caminho base onde as pastas "1", "2", ..., "10" estão localizadas
    output_path = "./traj_concatenate_aligned"  # Caminho da pasta de saída para a trajetória concatenada e alinhada
    concatenate_and_align_trajectories(base_path, output_path)
