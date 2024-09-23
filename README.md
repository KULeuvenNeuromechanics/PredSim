
*Use the table of contensts to easily navigate this README. Click on the three lines next to README.md just above this sentence.*

Predictive Simulations of Locomotion
============

<p align="center">
<img src="./Documentation/FiguresForDocumentation/PredSim%20logo.png" width="500" height="auto">
</p>

This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of locomotion. The simulation framework is [developed by Falisse *et al.*](#references). The implementation in this repo is aimed at letting you run simulations with your customized musculoskeletal models[*](#before-running-a-simulation).
If you want to reproduce specific published results, see [Replicate results](#replicate-results).

If you use PredSim for research, please cite 
> Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402.

> Lars D’Hondt, Antoine Falisse, Dhruv Gupta, Bram Van Den Bosch, Tom J. W. Buurke, Míriam Febrer-Nafría, Ines Vandekerckhove, Maarten Afschrift, and Friedl De Groote, PredSim: A Framework for Rapid Predictive Simulations of Locomotion, 10th IEEE RAS/EMBS International Conference for Biomedical Robotics and Biomechatronics (BioRob), Heidelberg, Germany, 2024

This repository is a work in progress.
1. ⚠️ **Troubleshooting Issues**:
   - If you encounter an error or get an unexpected output, let us know. Does the answer to your question/remark probably require code changes?
     - Yes -> use the [issues](https://github.com/KULeuvenNeuromechanics/PredSim/issues). Please look at the existing issues before creating a new one; you might not be the first person to have this problem.
     - No -> use the [discussions](https://github.com/KULeuvenNeuromechanics/PredSim/discussions). 
2. 🚧 **Suggesting Improvements**:
   - Feel free to suggest improvements. Submit a pull request with the changes, or post it in the discussion section on [Feature requests & Suggestions](https://github.com/KULeuvenNeuromechanics/PredSim/discussions/categories/feature-requests-suggestions).
3. 📯❔ **Sharing Best Practices and General Questions**:
   - To share best practices, publications, or more general questions, we encourage you to use the discussions section.



## Required software

To run this code you need to have the following softwares on your machine:

- MATLAB. [Statistics and Machine Learning Toolbox](https://nl.mathworks.com/products/statistics.html) and [Signal Processing Toolbox](https://nl.mathworks.com/products/signal.html) are required. [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html) is optional. The code has mainly been developed and tested on MATLAB 2021b, but is expected to run on any recent version.
- [OpenSim](https://simtk.org/projects/opensim) 4.3 or later. Older versions do not work.
- [CasADi](https://web.casadi.org/get/). The code has been tested on CasADi 3.5.5 and later.
- [Microsoft Visual Studio](https://visualstudio.microsoft.com/). In Visual Studio Installer, [select to include Desktop development with C++](./Documentation/FiguresForDocumentation/fig_MSVS.png). The code has been tested on MSVS Community 2015, 2017, 2019, and 2022.
- [CMake](https://cmake.org/download/). The code has been tested on CMake 3.22.0. Add CMake to your system Path.
- [Git](https://git-scm.com/download/win). The code has been tested on Git 2.40.0.windows.1. Add Git to your system Path.


## How to setup the code


1. Fork this repository to your github account. If you want the fork to be private, follow [these steps](./Documentation/PrivateForkPredSim.md) instead.
2. Clone the fork to your machine. If you have a computer with restricted permissions, make sure you have permission to run executables from the selected folder (For computers with KU Leuven BioMed Group policies, this is C:\GBW_MyPrograms\ ).
*Do not download the code as zip.*
3. Get the OpenSim API running on MATLAB. See [Setting up your Matlab Scripting Environment](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab#ScriptingwithMatlab-MatlabSetupSettingupyourMatlabScriptingEnvironment).
4. Add the location where you installed CasADi to the matlab search path (`addpath(genpath('C:/path/to/casadi'))`), or use S.solver.CasADi_path. 
5. Make sure the opensimAD submodule is installed. If PredSim\opensimAD\ is empty, open git command prompt, go to ...\PredSim\ , and run `git submodule update --init`.


After perfoming these steps, run the main script. (Expected run time is 40 minutes, depending on hardware.) If you don't receive any errors, your results should be the same as [the reference result](./Tests/ReferenceResults/Falisse_et_al_2022/Falisse_et_al_2022_paper.mat). If that is the case, you have succesfully intalled and set up the code. You are ready to do your own simulations.

## How to use the code

The code is written such that as a user you only have to interact with [*main.m*](./main.m).  

All user-defined [settings](./Documentation/SettingsOverview.md) are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings. 

The [examples folder](./Examples/) contains scripts illustrating how to use settings to:
- Predict a gait pattern resulting from a limited number of muscle synergies.
- Predict gait with assistance from an ankle exoskeleton.
- Run a sensitivity analysis.


### Before running a simulation

This code can automatically convert an OpenSim model to the external function used in the simulations. This greatly simplifies the process of going from a subject-specific model to a predictive simulation. Nevertheless, you should take care of the model you use since **not all OpenSim models are suported**: 
- Your model should not have locked joints. Locked joints would technically require having kinematic constraints, which is possible but makes the problem more complicated. Replace them with weld joints instead.
- Constraints on coordinates will be ignored (eg, coupling constraints).
- Using SimmSplines to describe coordinates (e.g. Yamaguchi knee model) is not supported as the implementation in OpenSim is not really compatible with algorithmic differentiation. Change them to Polynomials instead. GeometryPaths can contain SimmSplines. [_AdaptOpenSimModel.m_](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/AdaptOpenSimModel/AdaptOpenSimModel.m) takes care of changing present SimmSplines to polynomials.
- Your model needs to have contact elements that interact with the ground. Only *SmoothSphereHalfSpaceForce* contact forces are supported. You can use [_AdaptOpenSimModel.m_](./AdaptOpenSimModel/AdaptOpenSimModel.m) to add contact geometries and forces to your model. You can also scale the radius, stiffness and dissipation of the contact spheres.
- Your model can have any Hill-type muscle model, but it will be implemented as a [DeGroote-Fregly muscle](https://doi.org/10.1007/s10439-016-1591-9).
- Torque/force actuators of the class *ActivationCoordinateActuator* are supported. You can add actuators by running [_AdaptOpenSimModel.m_](./AdaptOpenSimModel/AdaptOpenSimModel.m). Actuators are not required.



### Required Settings

- **S.subject.name**: 
	- The name or code of the subject you are simulating. 
- **osim_path**: 
	- Path to the scaled opensim model of the subject.	
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

## References

1. Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim (MATLAB)

2. Falisse, A., Serrancolí, G., Dembia, C. L., Gillis, J., & De Groote, F. 2019 Algorithmic differentiation improves the computational efficiency of OpenSim-based trajectory optimization of human movement. PLOS ONE, 14(10), e0217730. https://doi.org/10.1371/journal.pone.0217730

3. Falisse A, Afschrift M, De Groote F. 2022 Modeling toes contributes to realistic stance knee mechanics in three-dimensional predictive simulations of walking. PLoS ONE 17(1): e0256311. https://doi.org/10.1371/journal.pone.0256311. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim_mtp (python)

4. D’Hondt, L., De Groote, F., & Afschrift, M. 2024 A dynamic foot model for predictive simulations of human gait reveals causal relations between foot structure and whole-body mechanics. PLOS Computational Biology, 20(6), e1012219. https://doi.org/10.1371/journal.pcbi.1012219. You can find the original repository here: https://github.com/Lars-DHondt-KUL/3dpredictsim/tree/four-segment_foot_model (MATLAB)
