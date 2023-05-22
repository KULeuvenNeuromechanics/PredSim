
*Use the table of contensts to easily navigate this README. Click on the three lines next to README.md just above this sentence.*

Predictive Simulations of Human Movement
============

This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of human gait. The simulation framework is [developed by Falisse *et al.*](#citations). The implementation in this repo is aimed at letting you run simulations with your customized musculoskeletal models[*](#before-running-a-simulations).
If you want to reproduce specific published results, you are adviced to visit the [repo corresponding to the paper](#citations).

This repository is a work in progress
- If you encounter an error or get an unexpected output, [check the issues](https://github.com/KULeuvenNeuromechanics/PredSim/issues). Please look at the existing issues before creating a new one, you might not be the first person to have this problem. Include the error messages in the issue.
- Feel free to suggest improvements. Submit a pull request with the changes, or create an issue with a desciption of the proposed change.

### Citations

Please cite the original authors if you use this simulation framework.

1. Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim (MATLAB)
 
2. Falisse A, Afschrift M, De Groote F. 2022 Modeling toes contributes to realistic stance knee mechanics in three-dimensional predictive simulations of walking. PLoS ONE 17(1): e0256311. https://doi.org/10.1371/journal.pone.0256311. You can find the original repository here: https://github.com/antoinefalisse/3dpredictsim_mtp (python)


### If you are familiar with a previous version of this simulation framework, here is a rundown of the changes

In general, hard coded variables were removed as much as possible and put in a settings structure the user can adapt. This results in:

- More easily use a different model (other number of joints or muscles)
- Make input-output a bit easier
- More easily change the (musculo-skeletal) parameters of the simulation

Besides that, the code was also adapted in accordance with recent (October 2020) adjustments in the original repository of Falisse et al.

- Adjusted version of the collocation scheme (now really being an orthogonal radau scheme)
- Support for parallel computing
- Formulation with opti

Lastly, seperate pieces of code were put together to streamline performing predictive simulations:

- Automatic conversion of an OpenSim model to the external function (executable called from the workflow)
- Perorming Muscle Analysis and polynomial fitting (integrated into the workflow)

## Required software

To run this code you need to have the following softwares on your machine:

- MATLAB. [Statistics and Machine Learning Toolbox](https://nl.mathworks.com/products/statistics.html) is required. [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html) is optional. The code has mainly been developed and tested on MATLAB 2021b, but is expected to run on any recent version.
- [OpenSim](https://simtk.org/projects/opensim) 4.3 or later. Older versions do not work.
- [CasADi](https://web.casadi.org/get/). The code has been tested on CasADi 3.5.5.
- [Microsoft Visual Studio](https://visualstudio.microsoft.com/). In Visual Studio Installer, [select to include Desktop development with C++](/FiguresForDocumentation/fig_MSVS.png). The code has been tested on MSVS Community 2015, 2017, 2019, and 2022.
- [CMake](https://cmake.org/download/). The code has been tested on CMake 3.22.0.


## How to setup the code


1. Fork this repository to your github account. If you want the fork to be private, follow [these steps](PrivateForkPredSim.md) instead.
2. Clone the fork to your machine. If you have a computer with restricted permissions, make sure you have permission to run executables from the selected folder (For computers with KU Leuven BioMed Group policies, this is C:\GBW_MyPrograms\ ).
*Do not download the code as zip.*
3. Get the OpenSim API running on MATLAB. See [Setting up your Matlab Scripting Environment](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab#ScriptingwithMatlab-MatlabSetupSettingupyourMatlabScriptingEnvironment).
4. In main.m, change [S.solver.CasADi_path](https://github.com/KULeuvenNeuromechanics/PredSim/blob/9fbbd43cf83617620e428d2c91f222c909a1349c/main.m#L84) to reflect the location where you installed CasADi. 
5. Make sure the opensimAD submodule is installed. If PredSim\opensimAD\ is empty, open git command prompt, go to ...\PredSim\ , and run `git submodule update --init`.


After perfoming these steps, run the main script. (Expected run time is 40 minutes, depending on hardware.) If you don't receive any errors, your results should be the same as https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Tests/Falisse_et_al_2022_Results. If that is the case, you have succesfully intalled and set up the code. You are ready to do your own simulations.

## How to use the code

The code is written such that as a user you only have to interact with [*main.m*](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/main.m).  

All user-defined settings are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings. 

### Before running a simulations

This code can automatically convert an OpenSim model to the external function used in the simulations. This greatly simplifies the process of going from a subject-specific model to a predictive simulation. Nevertheless, you should take care of the model you use since **not all OpenSim models are suported**: 
- Model should be 3D.
- Your model should not have locked joints. Locked joints would technically require having kinematic constraints, which is possible but makes the problem more complicated. Replace them with weld joints instead.
- Constraints on coordinates will be ignored (eg, coupling constraints).
- Using SimmSplines to describe coordinates (e.g. Yamaguchi knee model) is not supported as the implementation in OpenSim is not really compatible with algorithmic differentiation. Change them to Polynomials instead. GeometryPaths can contain SimmSplines.
- Your model needs to have contact elements that interact with the ground. Only *SmoothSphereHalfSpaceForce* contact forces are supported. You can use [_AdaptOpenSimModel.m_](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/AdaptOpenSimModel/AdaptOpenSimModel.m) to add contact geometries and forces to your model.
- Your model can have any Hill-type muscle model, but it will be implemented as a [DeGroote-Fregly muscle](https://doi.org/10.1007/s10439-016-1591-9).
- Torque/force actuators of the class *ActivationCoordinateActuator* are supported. You can add actuators by running [_AdaptOpenSimModel.m_](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/AdaptOpenSimModel/AdaptOpenSimModel.m). Actuators are not required.
- Ligament forces are not yet supported, but we plan to add them in the future.


### Required Settings

- **S.subject.name**: 
	- the name or code of the subject you are simulating.
- **osim_path**: 
	- path to the scaled opensim model of the subject.	
- **S.subject.save_folder**: 
	- path to the folder where you want to store the results of the OCP. If the folder does not exist yet on your machine, it will be created automatically.
- **S.subject.IG_selection**: 
	- either choose 'quasi-random' or give the path to a .mot file you want to use as initial guess.
- **S.subject.IG_selection_gaitCyclePercent**: 
	- if S.subject.IG_selection is a .mot file, S.subject.IG_selection_gaitCyclePercent is required. Here, specify what percent of gait cycle does the .mot file contain. For example, if the .mot file has 2 gait cycles, S.subject.IG_selection_gaitCyclePercent is 200.
- **S.solver.run_as_batch_job**: 
	- specify if the OCP is to be solved as a batch job (0: no, 1: yes). Batch processing requires the [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html).

### OptionalSettings

#### S.bounds

- **S.bounds.a.lower**: 
	- minimal muscle activation. Provide a number between 0 and 1. Default is *0.05* [double]
- **S.bounds.SLL.upper**: 
	- upper bound on left step length in meters. If not specified, no bound is implemented on left step length. 
- **S.bounds.SLR.upper**: 
	- upper bound on right step length in meters. If not specified, no bound is implemented on left step length.
- **S.bounds.dist_trav.lower**: 
	- lower bound on distance travelled over 1 gait cycle in meters. Note that if half gait cycle is being simulated, half of S.bounds.dist_trav.lower serves as the lower bound for total distance travelled. If not specified, no bound is implemented on left step length.
- **S.bounds.t_final.lower**: 
	- lower bound on final time in seconds. Default is *0.1* s [double].
- **S.bounds.t_final.upper**: 
	- upper bound on final time in seconds for full gait cycle simulation. Default is *2* s [double]. For half gait cycle simulation, half of this value gets implemented as upper bound for final time.
- **S.bounds.coordinates**: 
	- Cell array where 1st entry is dof name(s) , 2nd entry is its lower bound, and 3rd entry is its upper bound. Insert 'nan' or [] to lower bounds to only overwrite upper bounds, or vice versa. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',[],30} implements limit of -120° and 10° on knee angles, and default lower bound with 30° upper bound for pelvis_tilt. This setting changes the bounds of the optimization variables. When formulating the OCP, the variables are sclaed w.r.t. their bounds to improve conditioning. Changing these bounds can have a strong influence on convergence.
- **S.bounds.points**:
	- Cell array of structs where each cell defines a point. Points are used to define distanceconstraints. Each struct has the following fields:
		- body: name of a body in the OpenSim model [char].
		- point_in_body: xyz position of the point in the local frame of the body. Default is *[0, 0, 0]* [1x3 double].
		- name: name of the point. Default is name of the body [char]. 
- **S.bounds.distanceConstraints**:
	- Cell array of structs where each cell defines a constraint on the distance between two points. Each struct has the following fields:
		- point1: name of a point. If this is the name of a body in the OpenSim model, and no point with this name is defined, the origin of this body will be used [char]
		- point2: name of a point. If this is the name of a body in the OpenSim model, and no point with this name is defined, the origin of this body will be used [char]
		- direction: direction in which the distance  is constrained. Accepted inputs are: 1) any combination `x`, `y`, and `z`; and 2) `sagittal`, `coronal`, `frontal`, or `transverse` [char]. Note that for distances in one dimension (`point1 - point2`) the sign is kept 
		- lower_bound: lower bound on the distance, in m [double]. Default is no lower bound applied.
		- upper_bound: upper bound on the distance, in m [double]. Default is no upper bound applied.

#### S.metabolicE - metabolic energy

- **S.metabolicE.tanh_b**: 
	- hyperbolic tangeant smoothing factor used in the metabolic cost calculation. Default is *10* [double]
- **S.metabolicE.model**: 
	- the name of the metabolic energy model used. Default is [*Bhargava2004*](https://doi.org/10.1016/S0021-9290(03)00239-2) [char]. Currently only Bhargava2004 model has been implemented. Other options that will be added in the future are:
		- [*Umberger2003*](https://doi.org/10.1080/1025584031000091678)
		- [*Umberger2010*](https://doi.org/10.1098/rsif.2010.0084)
		- [*Uchida2016*](https://doi.org/10.1371/journal.pone.0150378)


#### S.misc - miscellanious

- **S.misc.gaitmotion_type**: 
	- type of gait simulation. Default is *HalfGaitCycle* [char]. Other option is *FullGaitCycle* [char]. Simulating a half gait cycle reduces computation time, but is limited to symmetric models. Post-processing will always reconstruct a full gait cycle starting at right heel strike.
- **S.misc.msk_geom_eq**: 
	- type of equation to approximate musculo-skeletal geometry (moment arm and muscle-tendon lengths wrt. joint angle). Default is *polynomials* [char]
- **S.misc.threshold_lMT_fit**: 
	- Threshold RMSE on muscle-tendon length to accept the polynomial fit. Default is *0.003* m [double]
- **S.misc.threshold_dM_fit**: 
	- Threshold RMSE on muscle-tendon moment arm to accept the polynomial fit. Default is *0.003* m [double]
- **S.misc.poly_order.lower**: 
	- minimal order of polynomial function. Default is *3* [double]
- **S.misc.poly_order.upper**: 
	- maximal order of polynomial function. Default is *9* [double]
- **S.misc.msk_geom_bounds**: 
	- Cell array where 1st entry is dof name(s) , 2nd entry is its lower bounds, and 3rd entry is its upper bounds. Insert nan to lower bounds to only overwrite upper bounds. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',nan,30} implements limit of -120 and 10 on knee angles, and default lower bund with 30 upper bound for pelvis_tilt. Defaults values are defined in the function [get_default_bounds_dummy_motion.m file](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/get_default_bounds_dummy_motion.m).
- **S.misc.visualize_bounds**: 
	- specify if bounds and initial guess are visualized (0 or 1). Default is *0* [double]
- **S.misc.dampingCoefficient**: 
	- damping coefficient of muscles. Default is *0.01* [double]. Used as damping value that is multiplied by the normalized muscle velocity, in the muscle velocity dependent term in calculation of normalized contractile element force of the muscle.
- **S.misc.constant_pennation_angle**: 
	- specify if pennation angle of the muscles is supposed to stay constant (0 or 1). Default is *0* [double]

#### S.post_process

- **S.post_process.make_plot**: 
	- boolean to plot post processing results (0 or 1). Default is *0*.
- **S.post_process.rerun**: 
	- boolean to rerun post-processing without solving OCP (0 or 1). Default is *0*. If this option is set to 1, one should specify the S.post_process.result_filename.
- **S.post_process.result_filename**: 
	- File name for results. Used for the name of .mat file that saves the results, diary of the OCP, and name of the .mot file of the output motion. When rerunning post-processing of an existing result, giving this file name is required. Default value is 
- **S.post_process.savename**: 
	- Type of savename to use if S.post_process.result_filename is empty. Defaults is *structured* [char]. This uses the name of the .mat file of results is used as <S.subject.name>_v\<n>. Where <S.subject.name> is defined in S.subject.name. n = 1 if <S.subject.name>_v1.mat does not exist. n is increased until n is found such that <S.subject.name>_v\<n>.mat does not exist. To change this structuring process, change its implementation in [run_pred_sim.m file](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/run_pred_sim.m). An alternative option is *datetime* [char], this uses <S.subject.name>\_\<yyyymmddTHHMMSS>. Where \<yyyymmddTHHMMSS> is the system date and time when creating the savename.

#### S.solver

- **S.solver.linear_solver**: 
	- solver algorithm used for the OCP. Default is *mumps* [char].
- **S.solver.tol_ipopt**: 
	- the power (10^-x) the tolerance ipopt has to reach before the OCP can be regarded as solved; a higher number gives a more precise answer. Default is *4* [double]
- **S.solver.max_iter**: 
	- maximal amount of iterations after wich the solver will stop. Default is *10000* [double]
- **S.solver.parallel_mode**: 
	- type of parallel computing. Default is *thread* [char].
- **S.solver.N_threads**: 
	- number of threads in parallel mode. Default is *4* [double]. When using batch computing, this value is overwritten with the number of threads assigned to each worker in you parallel cluster.
- **S.solver.N_meshes**: 
	- number of mesh intervals. Default is *50* [double] for S.misc.gaitmotion_type = HalfGaitCycle and *100* for FullGaitCycle

#### S.subject

- **S.subject.save_folder**: 
	- folder path to store the intermediate subject specific results (e.g. external function with skeletal dynamics, CasADi function with musculoskeletal geometry polynomials). If the folder does not exist, it is created automatically.
- **S.subject.mass**: 
	- mass of the subject in kilograms. Default is *[]* kilograms [double]. Default is empty, it will be overwritten by the mass extracted from the OpenSim model.
- **S.subject.IG_pelvis_y**: 
	- height from ground to pelvis, in meters. Default is *[]* m [double]. Default is empty, it will be overwritten by pelvis height extracted from the OpenSim model.
    - always used for the quasi-random initial guess
    - used for data-informed initial guess when `S.subject.adapt_IG_pelvis_y = 1;`
	- S.subject.IG_pelvis_y is also used to establish bounds on vertical pelvis position.
- **S.subject.adapt_IG_pelvis_y**: 
	- boolean to adjust the trajectory of height of pelvis from the ground for data-informed initial guess. Default is *0*. 0 means the trajectory will not be changed. If 1, the trajectory will be changed such that the average value of the trajectory is equal to s.subject.IG_pelvis_y.
- **S.subject.v_pelvis_x_trgt**: 
	- average velocity you want the model to have, in meters per second. Default is *1.25* m/s [double]
- **S.subject.IK_Bounds**: 
	- A .mot file that is used to define bounds on the kinematics. Default is *IK_Bounds_Default.mot* [char]. This file can be found in the OCP folder. This setting changes the bounds of the optimization variables. When formulating the OCP, the variables are sclaed w.r.t. their bounds to improve conditioning. Changing these bounds can have a strong influence on convergence.
- **S.subject.muscle_strength**: 
	- structure with [scaling factors for muscle strength](/FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). This scales the max muscle force of the active muscle force. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_strength = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_scale**: 
	- structure with [scaling factors for muscle passive stiffness](/FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_shift**: 
	- structure with [scaling factors for muscle passive stiffness shift](/FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). This property shifts the start of passive muscle force from normalized muscle length = 1 to the valuse specified in S.subject.muscle_pass_stiff_shift. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factore. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.tendon_stiff_scale**: 
	- structure with [scaling factors for tendon stiffnesses](/FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.mtp_type**: 
	- type of mtp joint. Default is *''* [char], which treats the mtp like any other joint. Select *'2022paper'* to use passive mtp joints whose kinematics do affect the crossing muscle-tendon units ([Falisse et al., 2022](#citations).
- **S.subject.scale_MT_params**: 
	- scale muscle tendon properties that are read from opensim model. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is what property you want to scale (FMo, lMo, lTs, alphao or vMmax), 3rd is the scale factor itself. If more than one scaling is to be performed, add 3 more inputs. For example, S.subject.scale_MT_params = {{'soleus_l','soleus_r'},'FMo',0.9,{'tib_ant_l'},'lTs',1.1} will scale max isometric force of both soleus by a factor of 0.9 and tendon slack length of tibialis anterior left by a scale of 1.1.
- **S.subject.damping_coefficient_all_dofs**: 
	- damping coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default is *0.1* Nms/rad [double]
- **S.subject.set_damping_coefficient_selected_dofs**: 
	- damping coefficient can be specified here for each coordinate individually. For example, S.subject.set_damping_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.12,{'knee_angle_l'},0.11} will put damping coefficient of both hip flexions to 0.12 and that of knee angle left to 0.11. If not defined here for a particular coordinateS.subject.damping_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.stiffness_coefficient_all_dofs**: 
	- stiffness coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default in *0* Nm/rad [double]
- **S.subject.set_stiffness_coefficient_selected_dofs**: 
	- stiffness coefficient can be specified here for each coordinate individually. For example, S.subject.set_stiffness_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.012,{'knee_angle_l'},0.011} will put stiffness coefficient of both hip flexions to 0.012 Nm/rad and that of knee angle left to 0.011 Nm/rad. If not defined here for a particular coordinate, S.subject.damping_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.set_limit_torque_coefficients_selected_dofs**: 
	- Default values of coordinate limit torques are defined in the function [get_default_coord_limit_torque_coefficients.m](https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/get_default_coord_limit_torque_coefficients.m). If values other than these are to be used, they can be specified here.
- **S.subject.base_joints_legs**:
	- Joint name that is the base of a leg, left and right. Default is 'hip' [char]
- **S.subject.base_joints_arms**:
	- Joint name that is the base of an arm, left and right. Default is 'acromial' [char]

#### S.weights

- **S.weights.E**: 
	- weight on metabolic energy rate. Default is *500* [double]
- **S.weights.E_exp**: 
	- exponent for the metabolic energy rate. Default is *2* [double]
- **S.weights.q_dotdot**: 
	- weight on joint accelerations. Default is *50000* [double]
- **S.weights.e_arm**: 
	- weight on arm excitations. Default is *10^6* [double]
- **S.weights.pass_torq**: 
	- weight on passive torques. Default is *1000* [double]
- **S.weights.pass_torq_includes_damping**: 
	- specify if damping torque = damping coefficient * coordinate velocity is to be included in the cost function (0 or 1). Default is 0 [double].
- **S.weights.a**: 
	- weight on muscle activations. Default is *2000* [double]
- **S.weights.slack_ctrl**: 
	- weight on slack controls. Default is *0.001* [double]

#### S.OpenSimADOptions
These settings are passed to OpenSimAD.

- **S.OpenSimADOptions.compiler**: 
	- command prompt argument for the compiler. [char]
	By default, PredSim will look for the most recent version that is installed in either `C:/Program Files/Microsoft Visual Studio/` or `C:/Program Files (x86)/Microsoft Visual Studio/`.
   	If you get an error about not finding a compiler, use this setting to specify your compiler:
       - Visual studio 2015: 'Visual Studio 14 2015 Win64'
       - Visual studio 2017: 'Visual Studio 15 2017 Win64'
       - Visual studio 2017: 'Visual Studio 16 2019'
       - Visual studio 2017: 'Visual Studio 17 2022'
- **S.OpenSimADOptions.verbose_mode**:
	- print outputs from windows command prompt to matlab command window (and log file). Default is *true* [bool].
- **S.OpenSimADOptions.verify_ID**:
	- verify the generated function versus the inverse dynamics tool in OpenSim. Default is *false* [bool].
- **S.OpenSimADOptions.jointsOrder**: 
	- If you want to choose the order of the joints outputs. Default is empty, which uses the joint order of the .osim file. [cell array of char]
- **S.OpenSimADOptions.coordinatesOrder**: 
	- If you want to choose the order of the coordinate outputs. Default is empty, which uses the coordinate order of the .osim file. [cell array of char]
	S.OpenSimADOptions.jointsOrder and S.OpenSimADOptions.coordinatesOrder are included in the settings to aid backward compatibility.
- **S.OpenSimADOptions.input3DBodyForces**:
	- add 3D force vectors that act on bodies. Default is empty. Needs further implementations before this can be used.
- **S.OpenSimADOptions.input3DBodyMoments**:
	- add 3D moment vectors that act on bodies. Default is empty. Needs further implementations before this can be used.
- **S.OpenSimADOptions.export3DPositions**:
	- export 3D position of points in bodies, in ground reference frame. Default is empty. Needs further implementations before this can be used.
- **S.OpenSimADOptions.export3DVelocities**:
	- export 3D velocity of points in bodies, in ground reference frame. Default is empty. Needs further implementations before this can be used.
- **S.OpenSimADOptions.exportGRFs**: 
	- Export total ground reaction forces of left and right side. Default is *true* [bool]
- **S.OpenSimADOptions.exportSeparateGRFs**: 
	- Export ground reaction forces of each contact element. Default is *true* [bool]
- **S.OpenSimADOptions.exportGRMs**: 
	- Export total ground reaction moments of left and right side. Default is *true* [bool]
- **S.OpenSimADOptions.exportContactPowers**: 
	- Export power due to vertical compression of each contact element. Default is *true* [bool]

