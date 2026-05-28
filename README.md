This repository contains the code needed to reproduce the results from [paper citation].  

# (1) Building and Compiling LAMMPS

My build of LAMMPS modifies the `bpm/spring` package to prevent failure in compression by modifying the `bond_bpm_spring.cpp` file in the BPM package. You can build LAMMPS from my GitHub [repository](https://github.com/lexiarlen/lammps.git). Alternatively, you can build the [latest version](https://github.com/lammps/lammps.git) of LAMMPS, replace the `bond_bpm_spring.cpp` file with the file from my repository, and recompile. This second option may be desirable because my repository is less likely to be up to date with the most recent version of LAMMPS.

See the LAMMPS [Git installation documentation](https://docs.lammps.org/Install_git.html) for directions on how to clone LAMMPS from Git and build it on your machine. You will need to configure LAMMPS with the bonded particle, molecule, and granular packages. The instructions below worked for me on my machine, but this may vary.

### 2. clone the GitHub repository

Either clone from the official LAMMPS repository or from my fork of LAMMPS. If you clone from the official LAMMPS repository, replace the `bond_bpm_spring.cpp` file with the file in my fork before building.

```bash
git clone -b release https://github.com/lammps/lammps.git mylammps
```


### 3. configure the build

Load additional modules

```bash
ml load openmpi/5.0.5; ml load cmake
```

From inside mylammps, run

```bash
cmake -S cmake -B build -D BUILD_MPI=on -D PKG_BPM=on -D PKG_MOLECULE=on -D PKG_GRANULAR=on 
```

If you previously configured LAMMPS and want to reconfigure the model (say with another package), remove the old build directory first:
```bash
rm -rf build
```

### 4. build lammps

```bash
make -C build -j $SLURM_CPUS_ON_NODE
```

### 5. add LAMMPS to your PATH

After compiling, the LAMMPS executable should be located at:

```bash
mylammps/build/lmp
````

You can run LAMMPS directly using the full path to the executable:

```bash
mpirun -np 1 /path/to/mylammps/build/lmp -in /path/to/input/script/in.lammps
```

Alternatively, you can add the build directory to your `PATH`:

```bash
export PATH=/path/to/mylammps/build:$PATH
```

Then LAMMPS can be run as:

```bash
mpirun -np 1 lmp -in /path/to/input/script/in.lammps
```

You can test the executable with:

```bash
lmp -h
```

Example LAMMPS input scripts can be found in:

```bash
mylammps/examples
```

# (2) Creating Packings

The analysis folder contains batch scripts to run packings for the narrow channel experiments and the size effect experiments. These scripts contain code to create a bidisperse packing of particles with the radius of the large particles $r_L$ and the radius of the small particles $r_S = 3/5 r_L$. The ratio of the number of large particles $N_L$ to the number of small particles $N_S$ chosen to be  $N_L/N_S=\frac{(1+\sqrt{5})}{4}$. The packing fraction is $\phi = 0.7$. 

After the LAMMPS scripts have been run, the batch script launches a postprocessing pipeline, running create_atomfiles.py. In the packing scripts, $r_L=1$, which speeds up computation. The particles are then rescaled by the radius_scaling parameter. Thus, if a packing of particles of size $r_L$ is desired in a domain with length $L$, then radius_scaling should be set to $r_L$ and the packing domain length $L_\text{pack} = L/\text{radius_scaling}$. Then, the rescaled data from the LAMMPS output file is written to an atomfile which LAMMPS can read in other scripts. 

# (3) Running Experiments
The code is organized into the two experiments, simple shear and narrow channel. Within each folder, there a subfolders for each experiment. In each subfolder, there is an analysis folder which contains the code to produce the data to make the figures in the paper. A brief outline of these folders is provided here:

Figures 3, A1: simple_shear/size_effect

Figure 4: simple_shear/kbe2sigma

Figure 5: simple_shear/kb2elastic

Figure 6: simple_shear/long_shear

Figures 7,8,9: narrow_channel/parameter_sweep

Figures C1, C2, C3, C4: narrow_channel/kc_checks


# (4) Analyzing experiments

All figures are generated using the Python notebook figures.ipynb in the main directory.

## Note on Sherlock SLURM paths

The included `.sbatch` scripts reflect the Sherlock HPC layout used for the simulations. In this setup, source code is stored under `/home/groups/earlew/arlenlex/P1-model-dev`, while large simulation outputs and intermediate data are stored under `/scratch/groups/earlew/arlenlex/P1-model-dev`. Users running on another system should edit or export `PROJECT_DIR`, `SCRATCH_PROJECT_DIR`, `MODULEFILES_DIR`, and `CONDA_SH` before submitting jobs.

