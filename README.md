Predictive simulations human movement
============

This repository is a fork of the code and data to generate three-dimensional muscle-driven predictive simulations of human gait as described in: Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. 

In this repository the original code of Falisse et al. was adapted to:

- Simulate walking with external support of an exoskeleton
- Simulate model with mtp joint
- Simulate Rajagopal model
- Make input-output a bit easier (in my opinion)

Besides that, the code was also adapted in accordance with recent (October 2020) adjustments in the original repository of Falisse et al.

- Adjusted version of the collocation scheme (now really being an orthogonal radau scheme).
- support for parallel computing
- formulation with opti

### Installation instructions

Clone this repo, including the opensim-AD submodule (clone with the `--recursive` flag, or run `git submodule update --init` later)

### Create all input for the simulations

When using a new/adapted musclulosketal model, you have to execute three steps to create the surrogate models and equations needed for optimization. An example of these steps are shown in the matlab script **./ConvertOsimModel/Example_PrepareOptimimzation.m** 

A summary of the steps:

#### 1. Polynomial fitting

The funciton FitPolynomials create a surrogate model, based on polynomial functions, to compute muscle-tendon lengths and moment arms from the joint kinematics. 

- First, we create a sample of joint angles (i.e. dummy motion) and run muscle analysis on this dummy motion to create a training dataset. Note that running the muscle analysis takes about 20 minutes. 
- Second, we fit the polynomials functions and save it in a spefici folder (input argument PolyFolder). You'll have to point to this folder using the settings *S.PolyFolder* when running the optimization

#### 2. Create casadi functions

In the next step we read the muscle-tendon parameters from the model, combine it with the polynomial functions and create a casadi function for most of the equations used in the optimization. This includes:

- Equations for metabolic energy
- Equations for muscle dynamics
- Equations for activation dynamaics
- Casadi version of the polynomial functions
- ....

These functions are saved in a specific folder (input argument CasadiFunc_Folders). You'll have to point to this folder using the settings *S.CasadiFunc_Folders* when running the optimization

#### 3. Automatically create .dll files

You have to provide a .cpp file that solves inverse dynamics with the current model you are using. Note that creating this .cpp file (with the correct modelling parameters) is still a manual step. The conversion from .cpp to .dll is automized in the function CreateDllFileFromCpp, which you can download here https://github.com/MaartenAfschrift/CreateDll_PredSim

#### 4. Run your simulation

You can now run your tracking or predictive simulations when pointing to the correct:

- Folder with polynomial functions: **S.PolyFolder**
- Casadi functions:  **S.CasadiFunc_Folders**
- .dll files including the file used:
  - the optimization: **S.ExternalFunction**
  - the post processing: **S.ExternalFunction2**



### Run Tracking and predictive simulations

You can run the tracking and predictive simulations using the functions in the folder **OCP**. This includes

**Gait 92 model**:  (https://simtk-confluence.stanford.edu/display/OpenSim/Gait+2392+and+2354+Models)

- f_PredSim_Gait92.m solves the predictive simulations with the gait 92 model
- f_TrackSim_Gait92.m solves the tracking simulations with the gait 92 model [not finished yet]
- f_LoadSim_Gait92.m post processing/analysis from the simulated states and controls (both tracking and predictive simulations). 

**Rajagopal model**:  (https://simtk.org/projects/full_body)

- f_PredSim_Rajagopal.m sovles the predictive simulations with the gait 92 model
- f_TrackSim_Rajagopal.m solves the tracking simulations with the Rajagopal model
- f_LoadSim_Rajagopal.m post processing/analysis from the simulated states and controls (both tracking and predictive simulations). 

This functions requires a matlabstructure (here S) with the settings for the optimization as input. The default settings for the optimization are added to this settings structure using the function *GetDefaultSettings(S)*. You can find an overview of the settings below.

Typically you well run the optimization with a specific setup structure and then analyse the simulation results. As an example:

```matlab
% settings....
S.ResultsFolder = 'NameFolderSimResults';
S.savename      = 'Resuls_DefaultGait92';

% her also other required settings

% Run simulation
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-process simulation results
```



### Plot output

You can use the function PlotResults_3DSim the create a default figure with a summary of the results. You can easily add multiple simulations to this figure. For example



```matlab

%.....
% simulation walking 1.25 m/s
S.v_tgt = 1.25;
S.savename = 'WalkingNormal'
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-process simulation results (saves results as S.savename with the extension _pp)

% simulate walking 0.5 m/s
S.v_tgt = 0.5;
S.savename = 'WalkingSlow'
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-process simulation results (saves results as S.savename with the extension _pp)

% plot figure to compare results (with the three optional input arguments here)
h = figure(); 	% new figure with handle
PlotResults_3DSim(fullfile(S.ResultsFolder,'WalkingNormal_pp.mat'),[1 0 0],'Normal',h,1.25,'speed'); 	% plot results of normal walking
PlotResults_3DSim(fullfile(S.ResultsFolder,'WalkingSlow_pp.mat'),[0 0 1],'Slow',h,0.5,'speed'); 	% plot results of slow walking on same figure

```



#### Settings- Required

- **PolyFolder**: Folder with the surrogate model for the muscle-tendon length and moment arms (from step 1 of the section "Create all input for the simulations"). This path to the folder is relative to the folder (./Polynomials) [string]
- **CasadiFunc_Folders**: Name of the folder with the casadifunctions (exported in step 2 of the section "Create all input for the simulations"). This path to the folder is relative to the folder (./CasadiFunctions) [string]
- **v_tgt**: imposed walking speed [double]
- **ModelName**: select type of musculoskeletal model. Currently the two options are (1) Gait92 or (2) Rajagopal [string] 
- **Mass**: mass of the subject in kg [double]
- **ExternalFunc:** Name of the .dll file used in the optimization (used for solving inverse dynamics). This file should be in the folder *./ExternalFunctions*. See step three of the section "Create all input for the simulations". [string].
- **ExternalFunc2:** Name of the .dll file used for post processing. This file should be in the folder *./ExternalFunctions*. See step three of the section "Create all input for the simulations" [string]
- **ResultsFolder**: folder the save the results [string]
- **Savename**: the of the results file [string]



#### Settings - optional

**Simulated motion**

- **Symmetric**: simulate symmetric motion (i.e. half a gait cycle), default is true [boolean]
- **Periodic**: simulate a periodic motion (i.e. full gait cycle), default is false [boolean]

**Settings formulation and solving NLP**

- **N**: number of mesh intervals (default is 50) [double]
- **NThreads**: number of threads for parallel computing (default is 2) [double]
- **linear_solver**: default is mumps [string]
- **tol_ipopt:** tolerance of ipopt solver
- **parallelMode**: default is thread

**Weights **

- **W.E**: weight metabolic energy rate (default is 500)
- **W.Ak**: weight joint accelerations (default is 50000)
- **W.ArmE**: weight arm excitations (default is 10^6)
- **W.passMom**: weight passive torques (default is 1000)
- **W.A**: weight muscle activations (default is 2000)
- **W.exp_E**: power metabolic energy (default is 2)
- **W.Mtp**: weight mtp excitations (default is 10^6)
- **W.u**: weight on excitations arms actuators (default is 0.001)
- **W.Lumbar: ** weight on miniizing lumbar activations (in Rajagopal model) (default is 10^5)

**Initial guess** (Note: I should improve this in the future)

- **IGmodeID**: initial guess based on (1)walking motion, (2) running motion, (3) previous solution in the *Results* folder, (4) previous solution in the *./IG/data folder* default is (1)

- **IGsel**: (1) quasi random initial guess (2) data-based initial guess (default is 2)

- **IKfile_guess**: relative path to IK file used for initial guess (used when IGsel = 2 and IGmodeID is 1 or 2). Default is *OpenSimModel\IK_Guess_Default.mat*

- **savename_ig**: name of the IK file used for initial guess (used when IGmodelID is 4). This file should be in *./IG/data folder*. [string]

- **ResultsF_ig:** name of Folder with IK file for setting *savename_ig* (see above) when IGmodelID is 3 [string].

- **IG_PelvisY**: height of the pelvis in the quasi-random initial guess (in m) [double]

**Adapting bounds**

- **IKfile_Bounds**: relative path to IK file used to determine bounds (i.e. 3 times ROM in IK file for all DOFs). Default is *OpenSimModel\IK_Guess_Default.mat*
- **Bounds.ActLower**: lower bound on all muscle activations
- **Bounds.ActLowerHip**: lower bound on activation of the hip muscles
- **Bounds.ActLowerKnee**: lower bound on activation of the knee muscles
- **S.Bounds.ActLowerAnkle**: lower bound on activation of the ankle muscles

**Kinematic constraints**

- **Constr.calcn:** minimal distance between calcneneus (origin) in the transversal plane. default is 0.09m [double]
- **Constr.toes:** minimal distance between toes (origin) in the transversal plane. default is 0.09m [double]
- **Constr.tibia:** minimal distance between tibia(?s) (origin) in the transversal plane. default is 0.09m [double]

**Exoskeleton control **

- **DataSet**: name of the folder with exoskeleton assistance profile (saved in the folder *./Data*) with a .mat file named *torque_profile.mat*. This mat file should contain the variables *time* and *torque* with the torque profile for one full stride.
- **ExoBool:** Boolean to select if you want to include the torque profile (i.e. use exoskeleton)
- **ExoScale:** scale factor for the torque profile.