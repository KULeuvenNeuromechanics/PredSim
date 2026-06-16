# Running PredSim on a VSC cluster

If your local machine has insufficient memory or if you need to run many
simulations, you can consider moving your workload to the HPC clusters.
KU Leuven provides compute resources to researchers in the [High Performance Computing](https://icts.kuleuven.be/sc/onderzoeksgegevens/hpc)
service. The HPC clusters of KU Leuven are part of the [Vlaams Supercomputer Centrum](https://www.vscentrum.be/) (VSC).

## Getting started on KU Leuven's HPC
In order to get started on the HPC, it is highly recommended to follow the
[HPC Introduction](https://hpcleuven.github.io/HPC-intro/) and/or
[Linux Introduction](https://hpcleuven.github.io/Linux-intro/) courses. To get
an overview of upcoming runs of these courses, see the [VSC Training](https://www.vscentrum.be/vsctraining)
calendar. Extensive documentation on how to use the HPC clusters in general is
available at https://docs.vscentrum.be/

> **_NOTE:_** In order to use MATLAB on the KU Leuven HPC clusters, you need
to be a member of the `lli_matlab` group because it is licensed software. After
creating your VSC account, you can request membership of this group by visiting
https://account.vscentrum.be/django/group/new and typing `lli_matlab` in the
text box of the "Join group" paragraph.

To interact with the cluster, go to https://ondemand.hpc.kuleuven.be/ and click on Login Server Shell Access. 
Check to which credit accounts you have access by running the `sam-balance` command.

In order to get a copy of the PredSim repository on the cluster, the following
commands can be used:

```bash
cd $VSC_DATA
git clone https://github.com/KULeuvenNeuromechanics/PredSim
cd PredSim
git submodule init
git submodule update
```

If this ran succesfully, you should now see PredSim in the Data directory. To check the Data directory, click on Files -> Data directory. 
<img width="920" height="56" alt="image" src="https://github.com/user-attachments/assets/7353381d-c5a6-4f2b-a3c4-d75f5d34879f" />

Make sure to checkout a branch that is suited to run on the Linux environment
of the cluster, for instance check that the `opensimAD_linux` submodule is
present.

## Running simulations on the HPC
### Step 1: edit the job script
Typically simulations are run in batch mode on the cluster, meaning you prepare
a so-called job script including all commands that need to be executed to run
a simulation. Below is an example of such a job script, which is also included
as the [run_simulation.slurm](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/Examples/VSC/run_simulation.slurm) file in the repository. For debugging or test
runs you can also execute the commands from the job script in a terminal on
the cluster.

```bash
#!/bin/bash -l

#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --account=<credit_account>
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=6000M

# Instead of installing dependencies yourself, you can simply load modules on
# the cluster. The versions are just examples, you can check for available
# modules with a command like "module spider OpenSim"
module load OpenSim/4.3-foss-2024a
module load OpenSimAD/20230213-foss-2024a
module load CMake/3.29.3-GCCcore-13.3.0
module load CasADi/3.7.0-gfbf-2024a

# See https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=7165&p=36596
export LD_PRELOAD="${EBROOTGCCCORE}/lib64/libstdc++.so.6"

# Update the Java library path in order enable using the OpenSim Java bindings
export _JAVA_OPTIONS="${_JAVA_OPTIONS} -Djava.library.path=$EBROOTOPENSIM/sdk/lib"

# Let MATLAB pick up OpenBLAS from a module
export BLAS_VERSION=${EBROOTOPENBLAS}/lib64/libopenblas.so
export LAPACK_VERSION=${EBROOTOPENBLAS}/lib64/libopenblas.so

# TODO Figure out if it useful to use more than a single thread
export OMP_NUM_THREADS=1
matlab -nodisplay -nosplash -singleCompThread -r "addpath('Examples/VSC'); 	run_on_VSC_cluster"
```

Replace the `<credit_account>` entry with your own.


> [!NOTE]
> Genius is currently largely unavailable due to maintenance. Use --cluster=wice 
> with --partition=batch_icelake instead. OpenSim and OpenSimAD modules are 
> available on wice. 

### Step 2: run the job script
Navigate to your [OnDemand Dashboard](https://ondemand.hpc.kuleuven.be/) > Login Server Shell Access > start a simulation
using the predefined settings:

```bash
cd $VSC_DATA/PredSim
sbatch run_simulation.slurm
```
To see the status of your job, execute `squeue -M ALL`, terminal output will
be written to the job output file (by default looking like `slurm-<jobid>.out`.

> [!NOTE]
> If you upload files from a Windows machine to the Linux cluster, you might
> receive errors saying your file contains DOS line breaks instead of UNIX
> line breaks. You can rectify this by running `dos2unix <fn>` on the cluster.
>

### A few useful Git commands
```squeue --clusters=_NAMECLUSTER_ --me``` Check queue (replace _Name Cluster_)

```git rev-parse --short HEAD``` Check current commit
Update to new commit (stash and keep changes):
```
git stash
git checkout master
git pull origin master
git submodule update --recursive
git stash pop
```
Updating to new commit without restoring changes
```
git restore main.m run_simulation.slurm
git checkout master
git pull origin master
git submodule update --recursive
```

### Troubleshooting
If you encounter an error like "Please add CasADi to the matlab search path", 
make sure you are using the latest version of run_on_VSC_cluster.m from the 
repository, which handles CasADi and OpenSim paths automatically. If the issue 
persists, you can find the correct CasADi path by running:

```
module load CasADi/3.7.0-gfbf-2024a
echo $EBROOTCASADI
```



