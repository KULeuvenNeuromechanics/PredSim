

## Running PredSim on a VSC cluster

KU Leuven provides compute resources to researchers in the [High Performance Computing](https://icts.kuleuven.be/sc/onderzoeksgegevens/hpc)
service. The HPC clusters of KU Leuven are part of the [Vlaams Supercomputer Centrum](https://www.vscentrum.be/) (VSC).
If your local machine has insufficient memory or if you need to run many
simulations, you can consider moving your workload to the HPC clusters.
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

In order to get a copy of the PredSim repository on the cluster, the following
commands can be used:

```bash
cd $VSC_DATA
git clone https://github.com/KULeuvenNeuromechanics/PredSim
cd PredSim
git submodule init
git submodule update
```

Make sure to checkout a branch that is suited to run on the Linux environment
of the cluster, for instance check that the `opensimAD_linux` submodule is
present.

Typically simulations are run in batch mode on the cluster, meaning you prepare
a so-called job script including all commands that need to be executed to run
a simulation. Below is an example of such a job script, which is also included
as the `run_simulation.slurm` file in the repository. For debugging or test
runs you can also execute the commands from the job script in a terminal on
the cluster.

```bash
#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=18
#SBATCH --nodes=1
#SBATCH --account=<credit_account>
#SBATCH --time=04:00:00

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
matlab -nodisplay -nosplash -singleCompThread -r "main"
```

Replace the `<credit_account>` entry with your own (you can check to which
credit accounts you have access by running the `sam-balance` command) and
submit the job script from your PredSim directory with `sbatch run_simulation.slurm`.
To see the status of your job, execute `squeue -M ALL`, terminal output will
be written to the job output file (by default looking like `slurm-<jobid>.out`).

> [!NOTE]
> If you upload files from a Windows machine to the Linux cluster, you might
> receive errors saying your file contains DOS line breaks instead of UNIX
> line breaks. You can rectify this by running `dos2unix <fn>` on the cluster.
