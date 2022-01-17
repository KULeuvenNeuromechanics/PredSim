Predictive simulations human movement
============

This repository is a fork of the code and data to generate three-dimensional muscle-driven predictive simulations of human gait as described in: Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. 

In this repository the original code of Falisse et al. was adapted to:

- More easily use a diffent model
- Make input-output a bit easier
- Change parameters of your musculo-skeletal model more easily 
- Simulate model with mtp joint
- Simulate Rajagopal model

Besides that, the code was also adapted in accordance with recent (October 2020) adjustments in the original repository of Falisse et al.

- Adjusted version of the collocation scheme (now really being an orthogonal radau scheme).
- support for parallel computing
- formulation with opti

### Code structure

- main
	- getDefaultSettings
	- preprocessing
		- osim2dll
		- get_model_info
		- read_and_scale_MTparameters
			- getMTparameters
			- scale_MTparameters
		- get_musculoskeletal_geometry 
			- muscle_analysis
			- polynomial_fit
		- update_model_info
	- createCasadiFunctions
	- OCP_formulations
	- post-processing


### Settings

All user-defined settings are stored in structure S. In main.m you have to specify the required settings and are free to change/add the optional settings. 

#### Required

- **S.subject.save_folder**: path to the folder where you want to store the results. If the path does not exist yet on your machine, it will be created automatically.
- **S.subject.name**: the name or code of the subject you are simulating.
- **S.subject.IG_selection**: either choose "quasi-random" or give the path to a .mot file you want to use as initial guess.
- **S.subject.IG_bounds**: give the path to a .mot file on which IG_bounds will be based.

#### Optional

**bounds**

- **S.bounds.a.lower**: minimal muscle activation. Provide a number between 0 and 1. Default is 0 [double]

- **S.bounds.calcn_dist.lower**: minimal distance between calcanei (origin) in the transversal plane. Default is 0.09 m [double]
- **S.bounds.toes_dist.lower**: minimal distance between toes (origin) in the transversal plane. Default is 0.10 m [double]
- **S.bounds.tibia_dist.lower**: minimal distance between tibiae (origin) in the transversal plane. Default is 0.11 m [double]


end
=========




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