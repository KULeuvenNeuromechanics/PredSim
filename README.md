### Important

This is a work in progress

- If you encounter an error or get an unexpected output, [check the issues](https://github.com/KULeuvenNeuromechanics/PredSim/issues). Please look at the existing issues before creating a new one, you might not be the first person to have this problem. Include error messages.
- Feel free to suggest improvements. Submit a pull request with the changes, or create an issue with desciption of the change.
- This readme is not up to date, but comments inside the code are.

*Use the table of contensts to easily navigate this README. Click on the three lines next to README.md just above this sentence.*

Predictive Simulations of Human Movement
============

This repository is a rework of the code and data to generate three-dimensional muscle-driven predictive simulations of human gait as described in: Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim

In general, hard coded variables were removed as much as possible and put in a settings structure the user can adapt. This results in:

- More easily use a diffent model (other number of joints or muscles)
- Make input-output a bit easier
- More easily change the (musculo-skeletal) parameters of the simulation

Besides that, the code was also adapted in accordance with recent (October 2020) adjustments in the original repository of Falisse et al.

- Adjusted version of the collocation scheme (now really being an orthogonal radau scheme)
- Support for parallel computing
- Formulation with opti

Lastly, seperate pieces of code where put together to streamline performing predictive simulations:

- Automatic conversion of an OpenSim model to the external function (executable called from the workflow)
- Perorming Muscle Analysis and polynomial fitting (integrated into the workflow)


## Code Structure

- main
	- getDefaultSettings
	- pre-processing
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

## Installation Instruction

## How to use the software

## Settings

All user-defined settings are stored in structure *S*. In *main.m* you have to specify the required settings and are free to change/add the optional settings. 

### Required

- **S.subject.save_results**: path to the folder where you want to store the results of the OCP. If the path does not exist yet on your machine, it will be created automatically.
- **S.subject.name**: the name or code of the subject you are simulating.
- **S.subject.IG_selection**: either choose "quasi-random" or give the path to a .mot file you want to use as initial guess.
- **S.subject.IG_bounds**: give the path to a .mot file on which IG_bounds will be based.

### Optional

#### S.bounds

- **S.bounds.a.lower**: minimal muscle activation. Provide a number between 0 and 1. Default is *0* [double]
- **S.bounds.calcn_dist.lower**: minimal distance between calcanei (origin) in the transversal plane. Default is *0.09* m [double]
- **S.bounds.toes_dist.lower**: minimal distance between toes (origin) in the transversal plane. Default is *0.10* m [double]
- **S.bounds.tibia_dist.lower**: minimal distance between tibiae (origin) in the transversal plane. Default is *0.11* m [double]
- **S.bounds.SLL.upper**: upper bound on left step length in meters. Default is *[]* m [double]
- **S.bounds.SLR.upper**: upper bound on right step length in meters. Default is *[]* m [double]
- **S.bounds.dist_trav.lower**: lower bound on distance travelled in meters. Default is *[]* m [double]
- **S.bounds.t_final.upper**: upper bound on final time in seconds. Default is *[]* s [double]
- **S.bounds.t_final.lower**: lower bound on final time in seconds. Default is *[]* s [double]

#### S.metabolicE - metabolic energy

- **S.metabolicE.tanh_b**: hyperbolic tangeant smoothing factor used in the metabolic cost calculation. Default is *100* [double]
- **S.metabolicE.model**: the name of the metabolic energy model used. Default is *Bhargava2004* https://doi.org/10.1016/S0021-9290(03)00239-2 [char]. Other options are:
	- *Umberger2003* https://doi.org/10.1080/1025584031000091678
	- *Umberger2010* https://doi.org/10.1098/rsif.2010.0084
	- *Uchida2016* https://doi.org/10.1371/journal.pone.0150378

#### S.misc - miscellanious

- **S.misc.v_max_s**: maximal contraction velocity identifier. Default is *0* [double] :warning: ***TO CHECK***
- **S.misc.gaitmotion_type**: type of gait simulation. Default is *HalfGaitCycle* [char] :warning: ***list all other options***
- **S.misc.msk_geom_eq**: type of equation to approximate musculo-skeletal geometry (moment arm and muscle-tendon lengths wrt. joint angle). Default is *polynomials* [char] :warning: ***list all other options***
- **S.misc.poly_order.lower**: minimal order of polynomial function. Default is *3* [double]
- **S.misc.poly_order.upper**: maximal order of polynomial function. Default is *9* [double]

#### S.post_process

- **S.post_process.make_plot**: boolean to plot post processing results. Default is *0* :warning: ***TO CHECK: possible to make it a string array?***
- **S.post_process.savename**: name used for saving the result files. Either choose your own naming or *structured*. Default is *structured* [char] :warning: ***add how the structured name looks like?***

#### S.solver

- **S.solver.linear_solver**: solver algorithm used for the OCP. Default is *mumps* [char] :warning: ***add different options***
- **S.solver.tol_ipopt**: the power (10^-x) the dual infeasibility has to reach before the OCP can be regarded as solved; a higher number gives a more precise answer. Default is *4* [double]
- **S.solver.max_iter**: maximal amount of itereations after wich the solver will stop. Default is *10000* [double]
- **S.solver.parallel_mode**: type of parallel computing. Default is *thread* [char]. Other options are:
	- ... :warning: ***add other option(s)***
- **S.solver.N_threads**: number of threads in parallel mode. Default is *4* [double]
- **S.solver.N_meshes**: number of mesh intervals. Default is *50* [double]

#### S.subject

- **S.subject.save_folder**: folder path to store the intermediate subject specific results (muscle analysis etc.). This setting is created automatically.
- **S.subject.mass**: mass of the subject in kilograms. Default is *[]* kilograms [double]. If left empty, it will be overwritten by the mass extracted from the OpenSim model.
- **s.subject.IG_pelvis_y**: height from the ground of the pelvis for the initial guess, in meters. Default is *[]* m [double]. I left empty, it will be overwritten by pelvis height extracted from the OpenSim model.
- **S.subject.v_pelvis_x_trgt**: average velocity you want the model to have, in meters per second. Default is *1.25* m/s [double]
- **S.subject.muscle_strength**: structure with scaling factors for muscle strength. Default is *[]*. 
	- :warning: ***ad something about the possibility to choose either musclewise or jointwise strength scaling***
	- :warning: ***retink how to write down the information about this setting***
- **S.subject.muscle_stiff**: structure with scaling factors for muscle stiffness. Default is *[]*. 
	- :warning: ***ad something about the possibility to choose either musclewise or jointwise scaling***
	- :warning: ***retink how to write down the information about this setting***
- **S.subject.muscle_sym**: structure discribing muscle symetries. Default is *[]*.
	- :warning: ***retink how to write down the information about this setting***
- **S.subject.tendon_stiff**: structure with tendon stiffnesses. Default is *[]*. 
	- :warning: ***ad something about the possibility to choose either musclewise or jointwise changing of tendon stiffness***
	- :warning: ***retink how to write down the information about this setting***
- **S.subject.mtp_type**: type of mtp joint you want to use in the simulation. Default is *[]* :warning: ***ask Lars to further clarify this setting***
- **S.subject.MT_params**: muscle tendon properties. Default is *[]*
	- :warning: ***rethink this; should there be a reference to the muscle_strength, muscle_stiff and tendon_stiff?***
- **S.subject.spasticity**: muscle spasticity. Default is *[]*
	- :warning: ***follow up***
- **S.subject.muscle_coordination**: muscle coordination.
	- :warning: ***currently commented out; follow this up***
	
#### S.weights

- **S.weights.E**: weight on metabolic energy rate. Default is *500* [double]
- **S.weights.E_exp**: exponent for the metabolic energy rate. Default is *2* [double]
- **S.weights.q_dotdot**: weight on joint accelerations. Default is *50000* [double]
- **S.weights.e_arm**: weight on arm excitations. Default is *10^6* [double]
- **S.weights.pass_torq**: weight on passive torques. Default is *1000* [double]
- **S.weights.a**: weight on muscle activations. Default is *2000* [double]
- **S.weights.e_mtp**: weight on mtp excitation. Default is *10^6* [double]
- **S.weights.slack_ctrl**: weight on slack controls. Default is *0.001* [double]
