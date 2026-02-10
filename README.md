*Use the table of contents to easily navigate this README. Click on the three lines next to README.md just above this sentence.*

Predictive Simulations of Locomotion
============

<p align="center">
<img src="./Documentation/FiguresForDocumentation/PredSim%20logo.png" width="500" height="auto">
</p>

This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of locomotion. The simulation framework is [developed by Falisse *et al.*](#references). The implementation in this repo is aimed at letting you run simulations with your customized musculoskeletal models[*](#before-running-a-simulation).
If you want to reproduce specific published results, see [Replicate results](#replicate-results).

If you use PredSim for research, please cite 
> Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402.

> Lars D’Hondt, Antoine Falisse, Dhruv Gupta, Bram Van Den Bosch, Tom J. W. Buurke, Míriam Febrer-Nafría, Ines Vandekerckhove, Maarten Afschrift, and Friedl De Groote, PredSim: A Framework for Rapid Predictive Simulations of Locomotion, 10th IEEE RAS/EMBS International Conference for Biomedical Robotics and Biomechatronics (BioRob), Heidelberg, Germany, 2024, pp. 1208-1213, [doi: 10.1109/BioRob60516.2024.10719735](http://dx.doi.org/10.1109/BioRob60516.2024.10719735).


This repository is a work in progress.
1. ⚠️ **Troubleshooting Issues**:
   - If you encounter an error or get an unexpected output, let us know. Does the answer to your question/remark probably require code changes?
     - Yes -> use the [issues](https://github.com/KULeuvenNeuromechanics/PredSim/issues). Please look at the existing issues before creating a new one; you might not be the first person to have this problem.
     - No -> use the [discussions](https://github.com/KULeuvenNeuromechanics/PredSim/discussions). 
2. 🚧 **Suggesting Improvements**:
   - Feel free to suggest improvements. Submit a pull request with the changes ([See Guidelines for contributing](https://github.com/KULeuvenNeuromechanics/PredSim/wiki/Guidelines-for-contributing-code)), or post it in the discussion section on [Feature requests & Suggestions](https://github.com/KULeuvenNeuromechanics/PredSim/discussions/categories/feature-requests-suggestions).
3. 📯❔ **Sharing Best Practices and General Questions**:
   - To share best practices, publications, or more general questions, we encourage you to use the discussions section.



## Running PredSim on a local machine

### Required software

To run this code you need to have the following softwares on your machine:

- MATLAB. [Statistics and Machine Learning Toolbox](https://nl.mathworks.com/products/statistics.html) and [Signal Processing Toolbox](https://nl.mathworks.com/products/signal.html) are required. [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html) is optional. The code has mainly been developed and tested on MATLAB 2021b, but is expected to run on any recent version.
- [OpenSim](https://simtk.org/projects/opensim) 4.3 or later. Older versions do not work.
- [CasADi](https://web.casadi.org/get/). The code has been tested on CasADi 3.5.5 and later.
- [Microsoft Visual Studio](https://visualstudio.microsoft.com/). In Visual Studio Installer, [select to include Desktop development with C++](./Documentation/FiguresForDocumentation/fig_MSVS.png). The code has been tested on MSVS Community 2015, 2017, 2019, and 2022.
- [CMake](https://cmake.org/download/). The code has been tested on CMake 3.22.0. Add CMake to your system Path.
- [Git](https://git-scm.com/download/win). The code has been tested on Git 2.40.0.windows.1. Add Git to your system Path.


### How to setup the code


1. Fork this repository to your github account. If you want the fork to be private, follow [these steps](./Documentation/PrivateForkPredSim.md) instead.
2. Clone the fork to your machine. If you have a computer with restricted permissions, make sure you have permission to run executables from the selected folder (For computers with KU Leuven BioMed Group policies, this is C:\GBW_MyPrograms\ ).
*Do not download the code as zip.*
3. Get the OpenSim API running on MATLAB. See [Setting up your Matlab Scripting Environment](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab#ScriptingwithMatlab-MatlabSetupSettingupyourMatlabScriptingEnvironment).
4. Add the location where you installed CasADi to the matlab search path (`addpath(genpath('C:/path/to/casadi'))`), or use S.solver.CasADi_path. 
5. Make sure the opensimAD submodule is installed. If PredSim\opensimAD\ is empty, open git command prompt, go to ...\PredSim\ , and run `git submodule update --init`.


After perfoming these steps, run the main script. (Expected run time is 40 minutes, depending on hardware.) If you don't receive any errors, your results should be the same as [the reference result](./Tests/ReferenceResults/Falisse_et_al_2022/Falisse_et_al_2022_paper.mat). If that is the case, you have succesfully intalled and set up the code. You are ready to do your own simulations.


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
matlab -nodisplay -nosplash -singleCompThread -r "addpath('Examples'); 	run_on_VSC_cluster"
```

Replace the `<credit_account>` entry with your own (you can check to which
credit accounts you have access by running the `sam-balance` command) and
submit the job script from your PredSim directory with `sbatch run_simulation.slurm`.
To see the status of your job, execute `squeue -M ALL`, terminal output will
be written to the job output file (by default looking like `slurm-<jobid>.out`.

> [!NOTE]
> If you upload files from a Windows machine to the Linux cluster, you might
> receive errors saying your file contains DOS line breaks instead of UNIX
> line breaks. You can rectify this by running `dos2unix <fn>` on the cluster.

If you have already completed these setup steps, you can simply navigate to 
your [OnDemand Dashboard](https://ondemand.hpc.kuleuven.be/) > Login Server Shell Access > start a simulation
using the predefined settings:
```bash
cd $VSC_DATA/PredSim
sbatch run_simulation.slurm
```


## How to use the code

The code is written such that as a user you only have to interact with [*main.m*](./main.m).  

All user-defined [settings](./Documentation/SettingsOverview.md) are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings. 


The [examples folder](./Examples/) contains scripts illustrating how to use settings to:
- Predict a gait pattern resulting from a limited number of muscle synergies.
- Predict gait with assistance from an ankle exoskeleton.
- Run a sensitivity analysis.
- Run simulations in batch on the VSC cluster.


### Required Settings

- **osim_path**: 
	- Path to the scaled opensim model of the subject. This code can automatically convert an OpenSim model to the external function used in the simulations. This greatly simplifies the process of going from a subject-specific model to a predictive simulation. Nevertheless, you should take care of the model you use since [**not all OpenSim models are suported**](./Documentation/ModelPreparation.md). 
- **S.subject.name**: 
	- The name or code of the subject you are simulating. 
- **S.misc.save_folder**: 
	- Path to the folder where you want to store the simulation results. If the folder does not exist yet on your machine, it will be created automatically.
- **S.solver.IG_selection**: 
	- Either choose 'quasi-random' or give the path to a .mot file you want to use as initial guess. In a quasi-random initial guess, the model is translated forward at the imposed velocity while all other coordinates are kept constant (vertical position of floating base is S.subject.IG_pelvis_y, others are 0). 
- **S.solver.IG_selection_gaitCyclePercent**: 
	- If S.solver.IG_selection is a .mot file, S.solver.IG_selection_gaitCyclePercent is required. Here, specify what percent of gait cycle does the .mot file contain. For example, if the .mot file has 2 gait cycles, S.solver.IG_selection_gaitCyclePercent is 200.

A [full overview of settings](./Documentation/SettingsOverview.md) is available in the documentation.

## Replicate results

Results of previous publications can be replicated by specifying one of the inputs for [initializeSettings](main.m#L18) and [S.subject.name](main.m#L23):
- `'Falisse_et_al_2022'`, see [reference 3](#references)
- `'DHondt_et_al_2024_3seg'`, see [reference 4](#references)
- `'DHondt_et_al_2024_4seg'`, see [reference 4](#references)

## Publications

For an overview of publications using PredSim, please check out the [publication list](./Documentation/PUBLICATIONS.md). 
Add your own references to highlight your work.

## References

1. Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim (MATLAB)

2. Falisse, A., Serrancolí, G., Dembia, C. L., Gillis, J., & De Groote, F. 2019 Algorithmic differentiation improves the computational efficiency of OpenSim-based trajectory optimization of human movement. PLOS ONE, 14(10), e0217730. https://doi.org/10.1371/journal.pone.0217730

3. Falisse A, Afschrift M, De Groote F. 2022 Modeling toes contributes to realistic stance knee mechanics in three-dimensional predictive simulations of walking. PLoS ONE 17(1): e0256311. https://doi.org/10.1371/journal.pone.0256311. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim_mtp (python)

4. D’Hondt, L., De Groote, F., & Afschrift, M. 2024 A dynamic foot model for predictive simulations of human gait reveals causal relations between foot structure and whole-body mechanics. PLOS Computational Biology, 20(6), e1012219. https://doi.org/10.1371/journal.pcbi.1012219. You can find the original repository here: https://github.com/Lars-DHondt-KUL/3dpredictsim/tree/four-segment_foot_model (MATLAB)
