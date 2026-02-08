#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import os

def check_overlaps(df):
    '''
    check for overlaps between atoms to avoid numerical blow ups
    '''
    coords = df[['x', 'y']].to_numpy()
    diameters = df["d"].to_numpy()
    n = len(coords)

    # compare each particle with every other particle (2d)
    all_overlaps = []
    for i in range(n):
        distances = np.sqrt(np.sum((coords[i] - coords)**2, axis=1))
        
        # get indices of particles that overlap where distance < diameter and not comparing with itself
        overlap_indices = np.where((distances < diameters[i]) & (distances > 0))[0]
        max_overlap = 0
        for j in overlap_indices:
            max_overlap = np.maximum(diameters[i] - distances[j], max_overlap)
        all_overlaps.append(max_overlap)
    abs_max_overlap = np.max(np.array(all_overlaps))

    if abs_max_overlap != 0:
        print(f"WARNING: overlapping particles with maximum overlap = {abs_max_overlap}")
    else:
        print("No overlaps found in dumped atoms.")
    return df

def main():
    if len(sys.argv) != 4:
        print("Usage: create_atomfiles.py path_to_dumped_atoms radius output_directory")
        sys.exit(1)

    _, path, radius, outputdir = sys.argv

    try:
        os.mkdir(outputdir)
    except FileExistsError:
        print('')
    except FileNotFoundError:
        print("")

    radius = float(radius)

    # read atom dump file
    df = pd.read_csv(path, sep=r'\s+', header=None, usecols = [0,1,2,3], names=["id", "type", "x", "y"], skiprows=9)
    df["d"] = (2.0-(df["type"].to_numpy()).astype(int))*0.4+0.6
    df.set_index("id")
    check_overlaps(df)
    num_atoms = len(df)
    vars = ["type", "x", "y"]

    # need to save the atom files for each variable in format readable for lammps 
    for var in vars:
        header_info = f"""# number of atoms below
{num_atoms}
# atom-id {var}
"""
        output_lines = [header_info]

        for index, row in df.iterrows():
            if var == "type":
                atom_line = f"{int(row["id"])} {int(row[var])}" 
            else:
                atom_line = f"{int(row["id"])} {row[var]*radius}" 
            output_lines.append(atom_line)

        output_content = "\n".join(output_lines)
        output_fname = os.path.join(outputdir, f'{var}.data')
        with open(output_fname, 'w') as file:
            file.write(output_content)
        print(f"File {output_fname} has been created successfully.")

if __name__ == "__main__":
    main()
