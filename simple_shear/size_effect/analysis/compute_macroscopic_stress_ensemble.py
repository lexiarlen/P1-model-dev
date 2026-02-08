from ovito.io import import_file
import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import os
import re
import matplotlib as mpl

def get_stress(filename, volume, phi, outputdir):
    pipeline = import_file(filename)

    total_sxy = []

    # get stress from each timestep with ovito
    for frame in range(pipeline.source.num_frames):
        data = pipeline.compute(frame)
        if 'c_peratom_stress[4]' in data.particles.keys():
            sxy = data.particles['c_peratom_stress[4]']
            total_sxy.append(np.sum(sxy))

    macroscopic_shear_stress = phi*np.array(total_sxy)/volume
    base = os.path.splitext(os.path.basename(filename))[0][:-7]
    np.save(os.path.join(outputdir, f'{base}_shear_array.npy'), macroscopic_shear_stress)
    return macroscopic_shear_stress

def main():
    # expecting 2 arguments: path_to_dumped_atoms output_directory
    if len(sys.argv) != 2:
        print("Usage: compute_macroscopic_stress.py path_to_dumped_atoms_dir")
        sys.exit(1)

    _, datadir = sys.argv

    outputdir = os.path.join(datadir, 'figures')
    os.makedirs(outputdir, exist_ok=True)

    ### hard coding file parsing ###
    Ls = np.array(['50', '100', '200', '300', '400', '500', '600'])
    seeds = np.array(['9834', '0825', '2893', '2283', '1923'])
    natoms_list = np.array(['873', '3496', '13987', '31473', '55952', '87426', '125893'])

    ### hard coding strain ###
    phi = 0.71

    shear_rate = 10**-7
    L_ini = 100e3
    volume = L_ini**2
    run_time = 1*3600 

    ### plotting ###

    fig, ax = plt.subplots()
    colors = ['indianred', 'sandybrown', 'gold', 'yellowgreen', 'skyblue', 'slateblue']

    for directory in sorted(os.listdir(datadir)):
        full_path = os.path.join(datadir, directory)
        if os.path.isdir(full_path) and directory.startswith("L"):
            k = 0
            Lval = (os.path.basename(directory))[1:]
            idx = (np.where(Ls == Lval)[0]).item()
            natoms = natoms_list[idx]
            macros = []
            print(f'working on directory = {directory} where L = {Lval} and natoms = {natoms}')
            for fname in sorted(glob.glob(os.path.join(os.path.join(datadir, directory), 'N*.lammps'))):
                print(f'individual file name is {fname}')
                macro = get_stress(fname, volume, phi, outputdir)
                macros.append(macro)
                k += 1
            macros = np.array(macros)
            mean_macro = np.mean(macros, axis = 0)
            var_macro = np.var(macros, axis = 0)
            std_error = np.sqrt(var_macro/k)
            time = np.linspace(0, run_time, len(macro))
            L = L_ini * (1+shear_rate*time)
            strain = (L-L_ini)/L_ini
            plt.plot(strain, mean_macro, color = colors[int(idx)], label = f'N = {natoms}', lw = 2)
            ax.fill_between(strain, mean_macro - std_error, mean_macro + std_error, color = colors[int(idx)], alpha = 0.3)
        plt.legend()
        plt.ylabel(r'$\tau$ [Pa]')
        plt.xlabel(r'$\gamma$')
        plt.ylim(-2500, 30000)
        plt.savefig(os.path.join(outputdir, 'stresses_w_dif_num_atoms.png'), dpi=300)


if __name__ == "__main__":
    main()
