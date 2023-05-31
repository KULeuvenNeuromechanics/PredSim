
*Use the table of contensts to easily navigate this README. Click on the three lines next to README.md just above this sentence.*

Predictive Simulations of Human Movement
============

![Predictive simulation of human walking](/FiguresForDocumentation/Falisse_et_al_2022_v1.gif)

This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of human gait. The simulation framework is [developed by Falisse *et al.*](#citations). The implementation in this repo is aimed at letting you run simulations with your customized musculoskeletal models[*](#before-running-a-simulations).
If you want to reproduce specific published results, you are adviced to visit the [repo corresponding to the paper](#citations).

This repository is a work in progress
- The main differences with the original framework developed by Falisse et. al. are listed on the wiki acompanying this repo, [here](https://github.com/KULeuvenNeuromechanics/PredSim/wiki/Current-repo-compared-to-Falisse-et.-al.-2022).
- If you encounter an error or get an unexpected output, [check the issues](https://github.com/KULeuvenNeuromechanics/PredSim/issues). Please look at the existing issues before creating a new one, you might not be the first person to have this problem. Include the error messages in the issue.
- Feel free to suggest improvements. Submit a pull request with the changes, or create an issue with a desciption of the proposed change.

### Citations

Please cite the original authors if you use this simulation framework.

1. Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim (MATLAB)
 
2. Falisse A, Afschrift M, De Groote F. 2022 Modeling toes contributes to realistic stance knee mechanics in three-dimensional predictive simulations of walking. PLoS ONE 17(1): e0256311. https://doi.org/10.1371/journal.pone.0256311. You can find the original repository here: https://github.com/antoinefalisse/predictsim_mtp (python)

## Required software

To run this code you need to have the following softwares on your machine:

- MATLAB. [Statistics and Machine Learning Toolbox](https://nl.mathworks.com/products/statistics.html) is required. [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html) is optional. The code has mainly been developed and tested on MATLAB 2021b, but is expected to run on any recent version.
- [OpenSim](https://simtk.org/projects/opensim) 4.3 or later. Older versions do not work.
- [CasADi](https://web.casadi.org/get/). The code has been tested on CasADi 3.5.5.
- [Microsoft Visual Studio](https://visualstudio.microsoft.com/). In Visual Studio Installer, [select to include Desktop development with C++](/FiguresForDocumentation/fig_MSVS.png). The code has been tested on MSVS Community 2015, 2017, 2019, and 2022.
- [CMake](https://cmake.org/download/). The code has been tested on CMake 3.22.0.


## How to setup the code


1. Fork this repository to your github account. If you want the fork to be private, follow [these steps](PrivateForkPredSim.md) instead.
2. Clone the fork to your machine. Make sure there are no spaces in the path. If you have a computer with restricted permissions, make sure you have permission to run executables from the selected folder (For computers with KU Leuven BioMed Group policies, this is C:\GBW_MyPrograms\ ).
*Do not download the code as zip.*
3. Get the OpenSim API running on MATLAB. See [Setting up your Matlab Scripting Environment](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab#ScriptingwithMatlab-MatlabSetupSettingupyourMatlabScriptingEnvironment).
4. In main.m, change [S.solver.CasADi_path](https://github.com/KULeuvenNeuromechanics/PredSim/blob/9fbbd43cf83617620e428d2c91f222c909a1349c/main.m#L84) to reflect the location where you installed CasADi. 
5. In main.m, change [S.Cpp2Dll.PathCpp2Dll_Exe](https://github.com/KULeuvenNeuromechanics/PredSim/blob/9fbbd43cf83617620e428d2c91f222c909a1349c/main.m#L115) to specify where you want to have the executable installed that will convert the OpenSim models to the external function. If you have a computer with KU Leuven GBW restrictions, be sure to have this path go into your 'C:\GBW_MyPrograms' folder.
> note: When running the code for the 1st time, it will download some dependencies. Mind your internet connection.
6. Make sure the opensimAD submodule is installed. If PredSim\opensimAD\ is empty, open git command prompt, go to ...\PredSim\ , and run `git submodule update --init`.


After perfoming these steps, run the main script. (Expected run time is 40 minutes, depending on hardware.) If you don't receive any errors, your results should be the same as https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Tests/Falisse_et_al_2022_Results. If that is the case, you have succesfully intalled and set up the code. You are ready to do your own simulations.

## How to use the code

The code is written such that as a user you only have to interact with [*main.m*](main.m).  

All user-defined settings are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings. 

### Before running a simulations

This code can automatically convert an OpenSim model to the external function used in the simulations. This greatly simplifies the process of going from a subject-specific model to a predictive simulation. Nevertheless, you should take care of the model you use since **not all OpenSim models are suported**: 
- Model should be 3D.
- Your model should not have locked joints. Locked joints would technically require having kinematic constraints, which is possible but makes the problem more complicated. Replace them with weld joints instead.
- Constraints on coordinates will be ignored (eg, coupling constraints).
- Using SimmSplines to describe coordinates (e.g. Yamaguchi knee model) is not supported as the implementation in OpenSim is not really compatible with algorithmic differentiation. Change them to Polynomials instead. GeometryPaths can contain SimmSplines. [_AdaptOpenSimModel.m_](AdaptOpenSimModel/AdaptOpenSimModel.m) takes care of changing present SimmSplines to polynomials.
- The kinematic chains starting at *acromial_l* and *acromial_r* will be interpreted as arms, legs start at *hip_l* and *hip_r*. A model is not required to have arms.
- Your model needs to have contact elements that interact with the ground. Only *SmoothSphereHalfSpaceForce* contact forces are supported. You can use [_AdaptOpenSimModel.m_](AdaptOpenSimModel/AdaptOpenSimModel.m) to add contact geometries and forces to your model.
- Your model can have any Hill-type muscle model, but it will be implemented as a [DeGroote-Fregly muscle](https://doi.org/10.1007/s10439-016-1591-9).
- Torque/force actuators of the class *ActivationCoordinateActuator* are supported. You can add actuators by running [_AdaptOpenSimModel.m_](AdaptOpenSimModel/AdaptOpenSimModel.m). Actuators are not required.
- Ligament forces are not yet supported, but we plan to add them in the future.
- If running simulation with different models of the same subject, be sure that the filename of the model is different for each model.

#### Required Settings

- **S.subject.name**: 
	- the name or code of the subject you are simulating.
- **S.subject.save_folder**: 
	- path to the folder where you want to store the results of the OCP. If the folder does not exist yet on your machine, it will be created automatically.
- **S.subject.IG_selection**: 
	- either choose 'quasi-random' or give the path to a .mot file you want to use as initial guess.
- **osim_path**: 
	- path to the scaled opensim model of the subject.
- **S.subject.IG_selection_gaitCyclePercent**: 
	- if S.subject.IG_selection is a .mot file, S.subject.IG_selection_gaitCyclePercent is required. Here, specify what percent of gait cycle does the .mot file contain. For example, if the .mot file has 2 gait cycles, S.subject.IG_selection_gaitCyclePercent is 200.
- **S.Cpp2Dll.PathCpp2Dll_Exe**:
	- path to folder with compiled [opensimAD](https://github.com/Lars-DHondt-KUL/opensimAD) files to create dll files from opensim model. It will be downloaded or updated automatically when required.
- **S.solver.run_as_batch_job**: 
	- specify if the OCP is to be solved as a batch job (0: no, 1: yes). Batch processing requires the [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html).

You can find the description of additional settings in [Settings.md](Settings.md)