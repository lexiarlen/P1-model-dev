# (1) Building and Compiling LAMMPS

Clone and install my fork of LAMMPS. This fork makes a simple change of preventing failure in compression and modifying the smoothing for the bond force in the ___.cpp file of the BPM package in LAMMPS. 

Go here [here](https://docs.lammps.org/Install_git.html) for directions on how to clone LAMMPS from git and build it on your machine. To build the version of LAMMPS with 

### 2. clone the github repository

```bash
git clone -b release https://github.com/lammps/lammps.git mylammps
```

### 3. configure the build

Load additional modules

```bash
ml load openmpi/5.0.5; ml load cmake
```

Configure the build

```bash
cmake -D BUILD_MPI=on -D PKG_BPM=on -D PKG_MOLECULE=on -D PKG_GRANULAR=on ../cmake
```

### 4. build lammps

```bash
make -j $SLURM_CPUS_ON_NODE
```

### 5. set up a symbolic link

```bash
mkdir -p /home/groups/earlew/arlenlex/lammps_lexi/v1/bin;
ln -sf /home/groups/earlew/arlenlex/mylammps/build/lmp \
       /home/groups/earlew/arlenlex/lammps_lexi/v1/bin/lmp
```

### 6. test lammps

follow the instructions below

# running

### Load the module:

```bash
ml purge
ml use /home/users/arlenlex/.local/modulefiles # questionable how accessible this is
ml load lammps_lexi/v1
```

### Then you can run your script using the executable lmp:

You can find example scripts to run in /home/groups/earlew/arlenlex/mylammps/examples.

# (2) Creating Packings

Creating packings. 

# (3) Running Experiments

# (4) Analyzing experiments

Need to create a conda environment with ovito. HUGE PAIN on hpc. pip install rather than conda. Virtual envs won't work. 