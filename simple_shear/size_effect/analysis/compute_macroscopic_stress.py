#!/usr/bin/env python3
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

    outputdir = os.path.join(datadir, 'output_files')
    os.makedirs(outputdir, exist_ok=True)

    # ----- hard coding shear calculation -----

    phi = 0.71
    L_ini = 100e3
    volume = L_ini**2

    # ------ iterate ------

    for directory in sorted(os.listdir(datadir)):
        full_path = os.path.join(datadir, directory)
        if os.path.isdir(full_path) and directory.startswith("L"):
            out_L_dir = os.path.join(outputdir, directory)
            os.makedirs(out_L_dir, exist_ok=True)

            for fname in sorted(glob.glob(os.path.join(full_path, 'N*.lammps'))):
                print(f'working on file {fname}', flush = True)
                macro = get_stress(fname, volume, phi)
                base = os.path.splitext(os.path.basename(fname))[0]

                outpath = os.path.join(out_L_dir, f'{base}_shear.npy')
                np.save(outpath, macro)

if __name__ == "__main__":
    main()
