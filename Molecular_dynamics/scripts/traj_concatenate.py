import os
import sys
import shutil
import warnings
from tqdm import tqdm  # Importar tqdm para a barra de progresso
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.DCD import DCDWriter

warnings.filterwarnings("ignore")


def concatenate_and_align_trajectories(base_path, output_path, name, number_folders):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    topology_path = os.path.join(base_path, f"1/{name}_1.prmtop")
    trajectory_paths = [os.path.join(base_path, f"{i}/{name}_{i}.dcd") for i in range(1,number_folders+1)]

    u = mda.Universe(topology_path, trajectory_paths[0])
    ref_atoms = u.select_atoms('backbone')

    ref_universe = mda.Universe(topology_path, trajectory_paths[0])
    ref_frame = ref_universe.trajectory[0]
    ref_positions = ref_universe.select_atoms('backbone').positions

    total_frames = sum([mda.Universe(topology_path, traj).trajectory.n_frames for traj in trajectory_paths])

    with DCDWriter(os.path.join(output_path, "all_traj_aligned.dcd"), n_atoms=u.atoms.n_atoms) as W:
        for traj_path in trajectory_paths:
            u.load_new(traj_path)
            for ts in tqdm(u.trajectory, total=u.trajectory.n_frames, desc=f"Processing {traj_path}"):
                align.alignto(u, ref_universe, select='backbone')
                W.write(u)

    shutil.copy(topology_path, os.path.join(output_path, f"{name}_1.prmtop"))

if __name__ == "__main__":
    name = sys.argv[1]
    number_folders = sys.argv[2]
    base_path = "./"
    output_path = "./traj_concatenate_aligned"
    concatenate_and_align_trajectories(base_path, output_path, name, int(number_folders))
