# obtain domain averaged shear stress from lammps dump files

from ovito.io import import_file
import numpy as np
import sys
import glob
import os

def get_stress(filename, volume, phi):
    pipeline = import_file(filename)

    total_sxy = []

    # get stress from each timestep with ovito
    for frame in range(pipeline.source.num_frames):
        data = pipeline.compute(frame)
        if 'c_peratom_stress[4]' in data.particles.keys():
            sxy = data.particles['c_peratom_stress[4]']
            total_sxy.append(np.sum(sxy))

    macroscopic_shear_stress = phi*np.array(total_sxy)/volume
    return macroscopic_shear_stress

def main():
    # expecting 2 arguments: path_to_dumped_atoms output_directory
    if len(sys.argv) != 2:
        print("Usage: compute_macroscopic_stress.py path_to_dumped_atoms_dir")
        sys.exit(1)

    _, datadir = sys.argv

    outputdir = os.path.join(datadir, 'figures')
    os.makedirs(outputdir, exist_ok=True)

    # ----- hard coding file parsing and shear calculation -----
    Ls = np.array(['100', '200', '300', '400', '500', '600'])
    natoms_list = np.array(['3496', '13987', '31473', '55952', '87426', '125893'])

    phi = 0.71
    L_ini = 100e3
    volume = L_ini**2

    # ------ iterate ------

    for directory in sorted(os.listdir(datadir)):
        full_path = os.path.join(datadir, directory)
        if os.path.isdir(full_path) and directory.startswith("L"):
            for fname in sorted(glob.glob(os.path.join(os.path.join(datadir, directory), 'N*.lammps'))):
                print(f'individual file name is {fname}')
                macro = get_stress(fname, volume, phi)
                base = os.path.splitext(os.path.basename(fname))[0][:-7]
                outpath = os.path.join(directory, f'{base}_shear.npy')
                np.save(outpath, macro)

if __name__ == "__main__":
    main()