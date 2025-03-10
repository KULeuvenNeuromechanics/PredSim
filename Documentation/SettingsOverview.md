
Quickly navigate to:
* [Required Settings](SettingsOverview.md#required-settings)
* [Optional Settings](SettingsOverview.md#optional-settings)
  * [S.bounds](SettingsOverview.md#sbounds)
  * [S.metabolicE - metabolic-energy](SettingsOverview.md#smetabolice---metabolic-energy)
  * [S.misc - miscellanious](SettingsOverview.md#smisc---miscellanious)
  * [S.post_process](SettingsOverview.md#spost_process)
  * [S.solver](SettingsOverview.md#ssolver)
  * [S.subject](SettingsOverview.md#ssubject)
  * [S.weights](SettingsOverview.md#sweights)
  * [S.OpenSimADOptions](SettingsOverview.md#sopensimadoptions)
  * [S.orthosis](SettingsOverview.md#sorthosis)

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



### Optional Settings

#### S.bounds

- **S.bounds.activation_all_muscles.lower**: 
    - minimal muscle activation. Provide a number between 0 and 1. Default is *0.05* [double]
- **S.bounds.activation_all_muscles.upper**: 
    - maximal muscle activation. Provide a number between 0 and 1. Default is *1* [double]
- **S.bounds.activation_selected_muscles**:
    - Cell array where 1st entry is muscle name(s) , 2nd entry is its lower bound, and 3rd entry is its upper bound. Insert 'nan' or [] to lower bounds to only overwrite upper bounds, or vice versa. For another bound, add 3 more entries.
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
- **S.bounds.points**:
    - Cell array of structs where each cell defines a point. Points are used to define distanceconstraints. Each struct has the following fields:
        - body: name of a body in the OpenSim model [char].
        - point_in_body: xyz position of the point in the local frame of the body. Default is *[0, 0, 0]* [1x3 double].
        - name: name of the point. Default is name of the body [char]. 
- **S.bounds.distanceConstraints**:
    - Cell array of structs where each cell defines a constraint on the distance between two points. Each struct has the following fields:
        - point1: name of a point. If this is the name of a body in the OpenSim model, and no point with this name is defined, the origin of this body will be used [char]
        - point2: name of a point. If this is the name of a body in the OpenSim model, and no point with this name is defined, the origin of this body will be used [char]
        - direction: direction in which the distance  is constrained. Accepted inputs are: 1) any combination `x`, `y`, and `z`; and 2) `sagittal`, `coronal`, `frontal`, or `transverse`. Default is *`xyz`* [char]. Note that for distances in one dimension (`point1 - point2`) the sign is kept.
        - lower_bound: lower bound on the distance, in m [double]. Default is no lower bound applied.
        - upper_bound: upper bound on the distance, in m [double]. Default is no upper bound applied.
- **S.bounds.Qs**: 
    - Cell array where 1st entry is dof name(s) , 2nd entry is its lower bound, and 3rd entry is its upper bound. In ° or m.
    Insert 'nan' or [] to lower bounds to only overwrite upper bounds, or vice versa. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',[],30} implements limit of -120° and 10° on knee angles, and default lower bound with 30° upper bound for pelvis_tilt. This setting changes the bounds of the optimization variables. When formulating the OCP, the variables are scaled w.r.t. their bounds to improve conditioning. Changing these bounds can have a strong influence on convergence.
- **S.bounds.Qdots**: 
    - Same as S.bounds.Qs, but for velocities.
- **S.bounds.Qdotdots**: 
    - Same as S.bounds.Qs, but for accelerations.
- **S.bounds.default_coordinate_bounds**:
    - Table with default values of bounds on Qs, Qdots, and Qdotdots. Default is *Default_Coordinate_Bounds.csv*. [string] Values in the file are assumed in rad or m.
- **S.bounds.Qdots_factor_RoM**:
    - Velocity bounds that are not given by another setting, will be taken symmetric and proportional to the range of motion. Default is *10* [double]
- **S.bounds.Qdotdots_factor_RoM**:
    - Acceleration bounds that are not given by another setting, will be taken symmetric and proportional to the range of motion. Default is *155* [double]
- **S.bounds.factor_IG_pelvis_ty.lower**:
    - Set lower bound op vertical position of floating base proportional to IG_pelvis_y. Default is *0.5* [double] Set to empty [] to not use this.
- **S.bounds.factor_IG_pelvis_ty.upper**:
    - Set upper bound op vertical position of floating base proportional to IG_pelvis_y. Default is *1.2* [double] Set to empty [] to not use this.

    Order of priority for coordinate bounds:
    1. Individual bounds from settings (S.bounds.Qs, S.bounds.Qdots, S.bounds.Qdotdots)
    2. Default bounds from table (S.bounds.default_coordinate_bounds)
    3. Read from model file. Qs: min and max coordinate values, Qdots: +/-10x coordinate range, Qdotdots: +/-155x coordinate range.


#### S.metabolicE - metabolic energy

- **S.metabolicE.tanh_b**: 
    - hyperbolic tangeant smoothing factor used in the metabolic cost calculation. Default is *10* [double]
- **S.metabolicE.model**: 
    - the name of the metabolic energy model used. Default is [*Bhargava2004*](https://doi.org/10.1016/S0021-9290(03)00239-2) [char]. Currently only Bhargava2004 model has been implemented. Other options that will be added in the future are:
        - [*Umberger2003*](https://doi.org/10.1080/1025584031000091678)
        - [*Umberger2010*](https://doi.org/10.1098/rsif.2010.0084)
        - [*Uchida2016*](https://doi.org/10.1371/journal.pone.0150378)


#### S.misc - miscellanious

- **S.misc.forward_velocity**:
    - Imposed forward velocity. Forward velocity is calculated as the average velocity in the coordinate that translates the floating base along the x-axis. Default is *1.25* [double].
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
- **S.misc.default_msk_geom_bound**:
    - file with default values for upper and lower bounds for approximating musculoskeletal geometry. Rotations are assumed in degrees, translations in meters. Default is *'default_msk_geom_bounds.csv'* [char].
    The provided file should be compatible with [`readtable`](https://mathworks.com/help/matlab/ref/readtable.html). The table should contain a column with coordinate names (header: name), a column with lower bounds, in degrees and meters (header: lower), and a column with upper bounds, in degrees and meters (header: upper). To set only the upper or lower bound of a coordinate, set the other one to `nan`;
- **S.misc.msk_geom_bounds**: 
    - Cell array where 1st entry is dof name(s) , 2nd entry is its lower bounds, and 3rd entry is its upper bounds. Insert nan to lower bounds to only overwrite upper bounds. For another bound, add 3 more entries. For example, {{'knee_angle_r','knee_angle_l'},-120,10,'lumbar_extension',nan,30} implements limit of -120° and 10° on knee angles, and default lower bound with 30° upper bound for lumbar_extension.
    > Order of priority for bounds:
    > 1. Individual bounds from settings (S.misc.msk_geom_bounds)
    > 2. Default bounds from table (S.misc.default_msk_geom_bound)
    > 3. Read from model file. Qs: min and max coordinate values
- **S.misc.msk_geom_n_samples**:
    - Number of samples for the dummy motion that is used to fit the approximated musculoskeletal geometry. Default is *5000* [double]

- **S.misc.visualize_bounds**: 
    - specify if bounds and initial guess are visualized (0 or 1). Default is *0* [double]
- **S.misc.dampingCoefficient**: 
    - damping coefficient of muscles. Default is *0.01* [double]. Used as damping value that is multiplied by the normalized muscle velocity, in the muscle velocity dependent term in calculation of normalized contractile element force of the muscle.
- **S.misc.constant_pennation_angle**: 
    - specify if pennation angle of the muscles is supposed to stay constant (0 or 1). Default is *0* [double]
- **S.misc.default_scaling_NLP**:
    - Filename with table that contains scale factors for Qs, Qdots, Qdotdots, and Moments in its columns. The first column should contain the coordinate names for its corresponding row. Default is *''*, i.e. scale factors are derived from bounds. [char]
- **S.misc.scaling_Qs**:
    - cell array of name-value pairs of coordinate names and the scale factor for the optimisation variables corresponding to their position.
- **S.misc.scaling_Qdots**:
    - cell array of name-value pairs of coordinate names and the scale factor for the optimisation variables corresponding to their velocity.
- **S.misc.scaling_Qdotdots**:
    - cell array of name-value pairs of coordinate names and the scale factor for the optimisation variables corresponding to their acceleration.
- **S.misc.scaling_moments**:
    - cell array of name-value pairs of coordinate names and the scale factor for the constraint violation on their moment equilibrium.
- **S.misc.git.local_hash**: 
    - hash of the local instance [char]. This is the identifier of the version of the code on your machine. You cannot change this setting.
- **S.misc.git.branch_name**: 
    - current branch of the local instance [char]. You cannot change this setting.
- **S.misc.git.remote_hash**: 
    - hash of the last commit on the remote [char]. This is the identifier of the latest version on the remote, i.e. GitHub. You cannot change this setting.
- **S.misc.computername**: 
    - name of the computer on which the simulation was run [char]. You cannot change this setting.
- **S.misc.save_folder**: 
    - path to folder to store the results. If the folder does not exist, it is created automatically. [char]
- **S.misc.result_filename**: 
    - File name for results. Used for the name of .mat file that saves the results, diary of the OCP, and name of the .mot file of the output motion. When rerunning post-processing of an existing result, giving this file name is required. Default value is empty.
- **S.misc.savename**: 
    - Type of savename to use if S.misc.result_filename is empty. Default is *structured* [char]. This sets S.misc.result_filename = <S.subject.name>_v\<n>. Where <S.subject.name> is defined in S.subject.name. n = 1 if <S.subject.name>_v1.mat does not exist. n is increased until n is found such that <S.subject.name>_v\<n>.mat does not exist. To change this structuring process, change its implementation in [run_pred_sim.m file](../run_pred_sim.m). An alternative option is *datetime* [char], this uses <S.subject.name>\_\<yyyymmddTHHMMSS>. Where \<yyyymmddTHHMMSS> is the system date and time when creating the savename.


#### S.post_process

- **S.post_process.rerun**: 
    - boolean to rerun post-processing without solving OCP (0 or 1). Default is *0*. If this option is set to 1, one should specify the S.misc.result_filename.
- **S.post_process.load_prev_opti_vars**:
    - load w_opt and reconstruct R before rerunning the post-processing. Advanced feature, for debugging only, you should not need this.

#### S.solver

- **S.solver.CasADi_path**:
    - Path to CasADi installation (top folder). By default, this will use the CasADi installation that is in the matlab search path.
- **S.solver.N_meshes**: 
    - number of mesh intervals. Default is *50* [double] for S.misc.gaitmotion_type = HalfGaitCycle and *100* for FullGaitCycle	
- **S.solver.run_as_batch_job**: 
    - specify if the OCP is to be solved as a batch job. Default is *false* [bool]. Batch processing requires the [Parallel Computing Toolbox](https://nl.mathworks.com/products/parallel-computing.html).
- **S.solver.batch_job_paths**:
    - if your simulation requires functions that are not inside the PredSim repo and you want to run this simulation as a batch job, then you have to include the path to the folder with these functions. Default is *{}* [cell array of char].
- **S.solver.par_cluster_name**:
    - name of the parallel cluster used to run batch jobs. Run `parallel.clusterProfiles` to see your available clusters. Default is *[]* [char] and will use your default cluster.
- **S.solver.parallel_mode**: 
    - type of parallel computing. Default is *thread* [char]. Other types are not supported.
- **S.solver.N_threads**: 
    - number of threads in parallel mode. Default is *4* [double]. When using batch computing, this value is overwritten with the number of threads assigned to each worker in your parallel cluster.
- **S.solver.nlpsol_options**:
    - Options to be passed to [CasADi's nonlinear program solver](https://web.casadi.org/api/html/d4/d89/group__nlpsol.html). Default is *[]* [struct].
- **S.solver.ipopt_options**:
    - Options to be passed to [ipopt](https://coin-or.github.io/Ipopt/OPTIONS.html#OPTIONS_REF). Default is *[]* [struct]. Overwrites S.solver.nlpsol_options.ipopt.
- **S.solver.linear_solver**: 
    - linear solver algorithm used by ipopt. Default is *mumps* [char]. Overwrites S.solver.ipopt_options.linear_solver.
- **S.solver.tol_ipopt**: 
    - the [ipopt convergence tolerance](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_tol) is set to 10^(-S.solver.tol_ipopt). A higher number gives a more precise answer, but takes more time. Default is *4* [double]. Overwrites S.solver.ipopt_options.tol.
- **S.solver.constr_viol_tol_ipopt**: 
    - the [ipopt constraint violation tolerance](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_constr_viol_tol) is set to 10^(-S.solver.constr_viol_tol_ipopt). A higher number gives a more precise answer, but takes more time. Default is *6* [double]. Overwrites S.solver.ipopt_options.constr_viol_tol.
- **S.solver.max_iter**: 
    - maximal amount of iterations after wich the solver will stop. Default is *10000* [double]. Overwrites S.solver.ipopt_options.max_iter.

#### S.subject

- **S.subject.mass**: 
    - mass of the subject in kilograms. Default is *[]* kilograms [double]. Default is empty, it will be overwritten by the mass extracted from the OpenSim model.
- **S.subject.IG_pelvis_y**: 
    - height from ground to pelvis, in meters. Default is *[]* m [double]. Default is empty, it will be overwritten by pelvis height extracted from the OpenSim model.
    - always used for the quasi-random initial guess
    - used for data-informed initial guess when `S.subject.adapt_IG_pelvis_y = 1;`
    - S.subject.IG_pelvis_y is also used to establish bounds on vertical pelvis position.
- **S.subject.adapt_IG_pelvis_y**: 
    - boolean to adjust the trajectory of height of pelvis from the ground for data-informed initial guess. Default is *0*. 0 means the trajectory will not be changed. If 1, the trajectory will be changed such that the average value of the trajectory is equal to S.subject.IG_pelvis_y.
- **S.subject.muscle_strength**: 
    - structure with [scaling factors for muscle strength](./FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). This scales the max muscle force of the active muscle force. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_strength = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_scale**: 
    - structure with [scaling factors for muscle passive stiffness](./FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.muscle_pass_stiff_shift**: 
    - structure with [scaling factors for muscle passive stiffness shift](./FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). This property shifts the start of passive muscle force from normalized muscle length = 1 to the valuse specified in S.subject.muscle_pass_stiff_shift. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factore. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.tendon_stiff_scale**: 
    - structure with [scaling factors for tendon stiffnesses](./FiguresForDocumentation/fig_muscle_tendon_properties_scaling.png). Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is the scale factor. If more than one scaling is to be performed, add 2 more inputs. For example, S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1} will scale both soleus by a factor of 0.9 and tibialis anterior left by a scale of 1.1.
- **S.subject.mtp_type**: 
    - type of mtp joint. Default is *''* [char], which treats the mtp like any other joint. Select *'2022paper'* to use passive mtp joints whose kinematics do affect the crossing muscle-tendon units ([Falisse et al., 2022](../README.md#references)).
- **S.subject.scale_MT_params**: 
    - scale muscle tendon properties that are read from opensim model. Default is *[]*, that is, no scaling. Input as a cell array where 1st input is the muscle(s) name, 2nd is what property you want to scale (FMo, lMo, lTs, alphao or vMmax), 3rd is the scale factor itself. If more than one scaling is to be performed, add 3 more inputs. For example, S.subject.scale_MT_params = {{'soleus_l','soleus_r'},'FMo',0.9,{'tib_ant_l'},'lTs',1.1} will scale max isometric force of both soleus by a factor of 0.9 and tendon slack length of tibialis anterior left by a scale of 1.1.
- **S.subject.damping_coefficient_all_dofs**: 
    - damping coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default is *0.1* Nms/rad [double]
- **S.subject.set_damping_coefficient_selected_dofs**: 
    - damping coefficient can be specified here for each coordinate individually. For example, S.subject.set_damping_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.12,{'knee_angle_l'},0.11} will put damping coefficient of both hip flexions to 0.12 and that of knee angle left to 0.11. If not defined here for a particular coordinateS.subject.damping_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.stiffness_coefficient_all_dofs**: 
    - stiffness coefficient for all coordinates (except coordinates connected to ground, generally pelvis (also called floating base)). Default in *0* Nm/rad [double]
- **S.subject.set_stiffness_coefficient_selected_dofs**: 
    - stiffness coefficient can be specified here for each coordinate individually. For example, S.subject.set_stiffness_coefficient_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},0.012,{'knee_angle_l'},0.011} will put stiffness coefficient of both hip flexions to 0.012 Nm/rad and that of knee angle left to 0.011 Nm/rad. If not defined here for a particular coordinate, S.subject.stiffness_coefficient_all_dofs will be used for that coordinate. Default is empty.
- **S.subject.set_stiffness_offset_selected_dofs**: 
    - position where moment of linear stiffness is zero can be specified here for each coordinate individually. For example, S.subject.set_stiffness_offset_selected_dofs = {{'hip_flexion_l','hip_flexion_r'},-0.3,} will offset the hip flexion stiffness such that their moment is zero at -0.3 rad hip flexion. Default is empty.
- **S.subject.default_coord_lim_torq_coeff**:
    - file with default coefficients for coordinate limit torques. Default is *'default_coord_lim_torq_coeff.csv'* [char].
    The provided file should be compatible with [`readtable`](https://mathworks.com/help/matlab/ref/readtable.html). The table should contain a column with coordinate names (header: name), 4 columns with stiffness coefficients, (headers: K_1, K_2, K_3, K_4), and 2 columns with offset coefficients (headers: theta_1, theta_2).
    Limit torques are calculated in function of coordinate value q as: `Tau = K(1)*exp(K(2)*(q-theta(2))) + K(3)*exp(K(4)*(q-theta(1)))`.
    Default coefficients are taken from *Anderson III, Frank Clayton. A dynamic optimization solution for a complete cycle of normal gait. The University of Texas at Austin, 1999.*
- **S.subject.scale_default_coord_lim_torq**:
    - scale factor for the amplitude of *all* default coordinate limit torques. Default is empty [double].
    All `K(1)` and `K(3)` values read from the file with default coefficients are multiplied by this factor.
    To scale the limit torques for individual degrees of freedom, use *S.subject.set_limit_torque_coefficients_selected_dofs*.
- **S.subject.set_limit_torque_coefficients_selected_dofs**: 
    - Set limit torque coefficients for a coordinate. Default is empty [cell array] with pattern {coordinate name(s) [char, cell array of chars], K [4x1 double], theta [2x1 double]}.
    For the specified coordinates, the coefficients from this setting take priority over the coefficients from the file with defaults (*S.subject.default_coord_lim_torq_coeff*).
    Note: Setting `K(1) = nan` will cause the limit torque for that coordinate to be excluded, even if it was given in the defaults.
- **S.subject.base_joints_legs**:
    - Joint name that is the base of a leg, left and right. Default is 'hip' [char]. Inputs of the form 'hip_r', {'hip_l'}, {'hip_r','hip_l'} are equivalent.
- **S.subject.base_joints_arms**:
    - Joint name that is the base of an arm, left and right. Default is 'acromial' [char]. Inputs of the form 'acromial_r', {'acromial_l'}, {'acromial_r','acromial_l'} are equivalent. Set to empty [] if the model does not have arms.
- **S.subject.stiffness_all_ligaments**:
    - Default stiffness model (i.e. force-length) used for ligaments. Default is [*ligamentGefen2002*](../ModelComponents/ligamentGefen2002.m) [char]
- **S.subject.set_stiffness_selected_ligaments**:
    - Use name-value pairs to use different stiffness models for specific ligaments. Default is *{'PlantarFascia',['plantarFasciaNatali2010'](../ModelComponents/plantarFasciaNatali2010.m)}* 
- **S.subject.synergies**:
	- boolean that indicates if muscle activations are controlled by synergies. Default is *0* (no synergies implemented).
	- When synergies are implemented, different variables need to be defined:
	- **S.subject.NSyn_r**: number of synergies for the right leg. This value needs to be defined.
	- **S.subject.NSyn_l**: number of synergies for the left leg. Default is equal to S.subject.NSyn_r. When simulating symmetric gait (i.e. S.misc.gaitmotion_type = 'HalfGaitCycle'), this is also set equal to S.subject.NSyn_r.
	- **S.subject.SynH_guess**: . Default is *0.1* [double]
	- **S.subject.SynW_guess**: . Default is *0.2* [double]
	- **S.subject.TrackSynW**: boolean that indicates if synergy weights are tracked. Default is *0* (no weights tracking).
	- **S.subject.TrackSynW_side**: indicates if weights are tracked for one or both legs. Possible options are: 'onlyLeft', 'onlyRight' or 'RightLeft'. Use this variable only if a full cycle is predicted (S.misc.gaitmotion_type = 'FullGaitCycle'). Default is *'RightLeft'* [char];
	- **S.subject.knownSynW_r**: synergy weights to be tracked are specified here, using the following form: {'muscleName1', weightArray1, {'muscleName2a', 'muscleName2b'}, weightArray2, etc.}. The 'weightArray' is an horizontal vector containing the muscle weight (between 0 and 1) on each synergy.
	- **S.subject.knownSynW_l**: same as 'S.subject.knownSynW_r' but for the left leg. Use this variable only if a full cycle is predicted (S.misc.gaitmotion_type = 'FullGaitCycle')
	- **S.subject.TrackSynW_NSyn_r**, **S.subject.TrackSynW_NSyn_l**: number of right and left synergy weights that are tracked for the muscles defined in 'S.subject.knownSynW_r' and 'S.subject.knownSynW_l'. This number may be different from 'S.subject.NSyn_r' or 'S.subject.NSyn_l'.

#### S.weights

- **S.weights.E**: 
    - weight on metabolic energy rate. Default is *500* [double]
- **S.weights.E_exp**: 
    - exponent for the metabolic energy rate. Default is *2* [double]
- **S.weights.q_dotdot**: 
    - weight on joint accelerations. Default is *50000* [double]
- **S.weights.e_torqAct**: 
    - weight on torque actuator excitations. Default is *10^6* [double]
- **S.weights.pass_torq**: 
    - weight on passive torques. Default is *1000* [double]. The passive torques term includes the coordinate limit torques, and potentially the coordinate damping torque.
- **S.weights.pass_torq_includes_damping**: 
    - specify if damping torque = damping coefficient * coordinate velocity is to be included in the cost function (0 or 1). Default is 0 [double].
- **S.weights.a**: 
    - weight on muscle activations. Default is *2000* [double]
- **S.weights.a_exp**: 
    - exponent for the muscle activation. Default is *2* [double]
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
       - Visual studio 2019: 'Visual Studio 16 2019'
       - Visual studio 2022: 'Visual Studio 17 2022'
- **S.OpenSimADOptions.verbose_mode**:
    - print outputs from windows command prompt to matlab command window (and log file). Default is *false* [bool].
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


#### S.orthosis
Interface to add custom orthoses to the simulation model.

- **S.orthosis.settings**:
    - Cell array of structs, where each struct contains the settings for one orthosis. Each orthosis should have a setting *function_name* [char] that refers to the function that describes the orthosis. See [functions in PredSim/WearableDevices for examples](../WearableDevices/). The other fields are passed to this function, so they can be used to set parameter values. See example functions and `help Orthosis` for how to define custom orthoses.

    

