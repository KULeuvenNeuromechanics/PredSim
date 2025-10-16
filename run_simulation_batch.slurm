#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --partition=batch
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
matlab -nodisplay -nosplash -singleCompThread -r "addpath('Examples'); 	run_on_VSC_cluster_batch"