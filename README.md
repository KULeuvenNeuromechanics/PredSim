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

- **S.subject.save_folder**: path to the folder where you want to store the results of the OCP. If the path does not exist yet on your machine, it will be created automatically.
- **S.subject.name**: the name or code of the subject you are simulating.
- **S.subject.IG_selection**: either choose "quasi-random" or give the path to a .mot file you want to use as initial guess.
- **S.subject.IG_selection_gaitCyclePercent**: If S.subject.IG_selection is a .mot file, S.subject.IG_selection_gaitCyclePercent is required. Here, specify what percent of gait cycle does the .mot file contain. For example, if the .mot file has 2 gait cycles, S.subject.IG_selection_gaitCyclePercent is 200.
- **S.solver.run_as_batch_job**: specify if the OCP is to be solved as a batch job (0 or 1). If solving as a batch, make sure batch processing is set up under parallel processing toolbox.
- **osim_path**: path to the scaled opensim model of the subject.

### Optional

#### S.bounds

- **S.bounds.a.lower**: minimal muscle activation. Provide a number between 0 and 1. Default is *0.05* [double]
- **S.bounds.calcn_dist.lower**: minimal distance between calcanei (origin) in the transversal plane. Default is *0.09* m [double]
- **S.bounds.toes_dist.lower**: minimal distance between toes (origin) in the transversal plane. Default is *0.10* m [double]
- **S.bounds.tibia_dist.lower**: minimal distance between tibiae (origin) in the transversal plane. Default is *0.11* m [double]
- **S.bounds.SLL.upper**: upper bound on left step length in meters. If not specified, no bound is implemented on left step length. 
- **S.bounds.SLR.upper**: upper bound on right step length in meters. If not specified, no bound is implemented on left step length.
- **S.bounds.dist_trav.lower**: lower bound on distance travelled over 1 gait cycle in meters. Note that if half gait cycle is being simulated, half of S.bounds.dist_trav.lower serves as the lower bound for total distance travelled. If not specified, no bound is implemented on left step length.
- **S.bounds.t_final.lower**: upper bound on final time in seconds. Default is *0.1* s[double].
- **S.bounds.t_final.upper**: lower bound on final time in seconds for full gait cycle simulation. Default is *2* s[double]. For half gait cycle simulation, half of this value gets implemented as upper bound for final time.
- **S.bounds.coordinates**: Cell array where 1st entry is dof name(s) , 2nd entry is its lower bounds, and 3rd entry is its upper bounds. Insert nan to lower bounds to only overwrite upper bounds. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',nan,30} implements limit of -120 and 10 on knee angles, and default lower bund with 30 upper bound for pelvis_tilt.

#### S.metabolicE - metabolic energy

- **S.metabolicE.tanh_b**: hyperbolic tangeant smoothing factor used in the metabolic cost calculation. Default is *100* [double]
- **S.metabolicE.model**: the name of the metabolic energy model used. Default is *Bhargava2004* https://doi.org/10.1016/S0021-9290(03)00239-2 [char]. Currently only Bhargava2004 model has been implemented. Other options are:
	- *Umberger2003* https://doi.org/10.1080/1025584031000091678
	- *Umberger2010* https://doi.org/10.1098/rsif.2010.0084
	- *Uchida2016* https://doi.org/10.1371/journal.pone.0150378


#### S.misc - miscellanious

- **S.misc.gaitmotion_type**: type of gait simulation. Default is *HalfGaitCycle* [char]. Other option is *FullGaitCycle* [char]
- **S.misc.msk_geom_eq**: type of equation to approximate musculo-skeletal geometry (moment arm and muscle-tendon lengths wrt. joint angle). Default is *polynomials* [char]
- **S.misc.threshold_lMT_fit**: Threshold RMSE on muscle-tendon length to accept the polynomial fit. Default is *0.003* m [double]
- **S.misc.threshold_dM_fit**: Threshold RMSE on muscle-tendon momnt arm to accept the polynomial fit. Default is *0.003* m [double]
- **S.misc.poly_order.lower**: minimal order of polynomial function. Default is *3* [double]
- **S.misc.poly_order.upper**: maximal order of polynomial function. Default is *9* [double]
- **S.misc.msk_geom_bounds**: Cell array where 1st entry is dof name(s) , 2nd entry is its lower bounds, and 3rd entry is its upper bounds. Insert nan to lower bounds to only overwrite upper bounds. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',nan,30} implements limit of -120 and 10 on knee angles, and default lower bund with 30 upper bound for pelvis_tilt. Defaults values are defined in the function get_default_bounds_dummy_motion.m file (https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/get_default_bounds_dummy_motion.m).
- **S.misc.visualize_bounds**: specify if bounds and initial guess are visualized (0 or 1). Default is *0.003* m [double]
- **S.misc.dampingCoefficient**: damping coefficient of muscles. Default is *0.01* [double]. Used as damping value that is multiplied by the normalized muscle velocity, in the muscle velocity dependent term in calculation of normalized contractile element force of the muscle.
- **S.misc.constant_pennation_angle**: specify if pennation angle of the muscles is supposed to stay constant (0 or 1). Default is *0* [double]

#### S.post_process

- **S.post_process.make_plot**: boolean to plot post processing results (0 or 1). Default is *0*.
- **S.post_process.rerun**: boolean to rerun post-processing without solving OCP. Default is *0*.
- **S.post_process.result_filename**: File name for results. Used for the name of .mat file that saves the results, diary of the OCP, and name of the .mot file of the output motion.
- **S.post_process.savename**: If S.post_process.result_filename is empty, S.post_process.savename is used. Defaults is *structured* [char]. This uses the name of the .mat file of results is used as <subject name>_v<n>. Where <subjenctname> is defined in S.subject.name. n = 1 if <subject name>_v1.mat does not exist. n is increased until n is found such that <subject name>_v<n>.mat does not exist. To change this structuring process, change its implementation in run_pred_sim.m file (https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/run_pred_sim.m).

#### S.solver

- **S.solver.linear_solver**: solver algorithm used for the OCP. Default is *mumps* [char]
- **S.solver.tol_ipopt**: the power (10^-x) the dual infeasibility has to reach before the OCP can be regarded as solved; a higher number gives a more precise answer. Default is *4* [double]
- **S.solver.max_iter**: maximal amount of itereations after wich the solver will stop. Default is *10000* [double]
- **S.solver.parallel_mode**: type of parallel computing. Default is *thread* [char].
- **S.solver.N_threads**: number of threads in parallel mode. Default is *4* [double]
- **S.solver.N_meshes**: number of mesh intervals. Default is *50* [double]

#### S.subject

- **S.subject.save_folder**: folder path to store the intermediate subject specific results (muscle analysis etc.). If the folder does not exist, it is created automatically.
- **S.subject.mass**: mass of the subject in kilograms. Default is *[]* kilograms [double]. Default is empty. It will be overwritten by the mass extracted from the OpenSim model.
- **s.subject.IG_pelvis_y**: height from the ground of the pelvis for the quasi-random initial guess, in meters. Default is *[]* m [double]. Default is empty, it will be overwritten by pelvis height extracted from the OpenSim model.
- **s.subject.adapt_IG_pelvis_y**: boolean to adjust the trajectory of height of pelvis from the ground for data-informed initial guess. Default is *0*. 0 means the trajectory will not be changed. If 1, the trajectory will be changed such that the average value of the trajectory is equal to s.subject.IG_pelvis_y.
- **S.subject.v_pelvis_x_trgt**: average velocity you want the model to have, in meters per second. Default is *IK_Bounds_Default.mot* [char]. This file can be found in the OCP folder.
- **S.subject.IK_Bounds**: A .mot file that is used to define bounds on the kinematics. Default is *1.25* m/s [double]
- **S.subject.muscle_strength**: structure with scaling factors for muscle strength. This scales the max muscle force of the active muscle force. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_strength = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_scale**: structure with scaling factors for muscle passive stiffness. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_shift**: structure with scaling factors for muscle passive stiffness shift. This property shifts the start of passive muscle force from normalized muscle length = 1 to the valuse specified in S.subject.muscle_pass_stiff_shift. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factore. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.tendon_stiff_scale**: structure with scaling factors for tendon stiffnesses. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.mtp_type**: type of mtp joint. Default is *''* [char]. The current implementation is the same as 2022 paper.
- **S.subject.scale_MT_params**: scale muscle tendon properties that are read from opensim model. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is what property you want to scale (FMo, lMo, lTs, alphao or vMmax), 3rd is the scale factor itself. If more than one scaling is to be performed, add 3 more inputs. For example, S.subject.scale_MT_params = {{'soleus_l','soleus_r'},'FMo',0.9,{'tib_ant_l'},'lTs',1.1} will scale max isometric force of both soleus by a factor of 0.9 and tendon slack length of tibialis anterior left by a scale of 1.1.
- **S.subject.damping_coefficient_all_dofs**: damping coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default is *0.1* Nms/rad [double]
- **S.subject.set_damping_coefficient_selected_dofs**: damping coefficient can be specified here for each coordinate individually. For example, S.subject.set_damping_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.12,{'knee_angle_l'},0.11} will put damping coefficient of both hip flexions to 0.12 and that of knee angle left to 0.11. If not defined here for a particular coordinateS.subject.damping_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.stiffness_coefficient_all_dofs**: stiffness coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default in *0* Nm/rad [double]
- **S.subject.set_stiffness_coefficient_selected_dofs**: stiffness coefficient can be specified here for each coordinate individually. For example, S.subject.set_stiffness_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.012,{'knee_angle_l'},0.011} will put stiffness coefficient of both hip flexions to 0.012 and that of knee angle left to 0.011. If not defined here for a particular coordiante, S.subject.damping_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.set_limit_torque_coefficients_selected_dofs**: Default values of coordinate limit torques are defined in the function get_default_coord_limit_torque_coefficients.m (https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/get_default_coord_limit_torque_coefficients.m). If values other than these are to be used, they can be specified here.
	
#### S.weights

- **S.weights.E**: weight on metabolic energy rate. Default is *500* [double]
- **S.weights.E_exp**: exponent for the metabolic energy rate. Default is *2* [double]
- **S.weights.q_dotdot**: weight on joint accelerations. Default is *50000* [double]
- **S.weights.e_arm**: weight on arm excitations. Default is *10^6* [double]
- **S.weights.pass_torq**: weight on passive torques. Default is *1000* [double]
- **S.weights.pass_torq_includes_damping**: specify if damping torque = damping coefficient * coordinate velocity is to be included in the cost function (0 or 1). Default is 1 [double].
- **S.weights.a**: weight on muscle activations. Default is *2000* [double]
- **S.weights.slack_ctrl**: weight on slack controls. Default is *0.001* [double]

#### S.Cpp2Dll - These settn=ings are not used during creation of the external function, and not during the OCP
- **S.Cpp2Dll.compiler**: select compiler for cpp projects. For example, 'Visual Studio 14 2015 Win64' or 'Visual Studio 15 2017 Win64'. Default is *Visual Studio 15 2017 Win64* [char]
- **S.Cpp2Dll.PathCpp2Dll_Exe**: Path with exectuables to create .cpp file. You can use the function S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe(ExeDir) to download this exectuable with the input 'ExeDir' to folder in which you want to install the executable. The output argument of this function gives you the path to the folder with the exectutable. Default is empty.
- **S.Cpp2Dll.export3DSegmentOrigins**: Export 3D segment origins. Default is S.Cpp2Dll.export3DSegmentOrigins = {'calcn_r', 'calcn_l', 'femur_r', 'femur_l', 'hand_r','hand_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l'};
- **S.Cpp2Dll.jointsOrder**: If you want to choose the order of the joints outputs. Default is empty, which uses the joint order of the .osim file.
- **S.Cpp2Dll.coordinatesOrder**: If you want to choose the order of the coordinate outputs. Default is empty, which uses the coordinate order of the .osim file.
- **S.Cpp2Dll.exportGRFs**: Export total GRFs (0 or 1). If True, right and left 3D GRFs (in this order) are exported. Set False or do not pass as argument to not export those variables. Default is 1.
- **S.Cpp2Dll.exportSeparateGRFs**: Export separate GRFs (0 or 1). If True, right and left 3D GRFs (in this order) are exported for each of the contact spheres. Set False or do not pass as argument to not export those variables. Default is 1.
- **S.Cpp2Dll.exportGRMs**: Export GRMs (0 or 1). If True, right and left 3D GRMs (in this order) are exported. Set False or do not pass as argument to not export those variables. Default is 1.
- **S.Cpp2Dll.exportContactPowers**: Export contact sphere vertical deformation power (0 or 1). If True, right and left vertical deformation power of all contact spheres are exported. Set False or do not pass as argument to not export those variables. Default is 1.
- **S.Cpp2Dll.verbose_mode**: Verbose mode (0 or 1). 0: only warnings and errors, 1: all information on building .dll file
