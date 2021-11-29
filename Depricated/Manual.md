# 3D predictive simulation



## Step-by-step guide

### Polynomial approximation

If you want to use a new model, you first have to approximate the muscle tendon length and moment arms as function of the joint kinematics using polynomials. We estimate the coefficient of these polynomials on a training dataset using the scripts in the folder Polynomials. For example for the model of the first subject of Katie Poggensee (i.e. exoskeleton torque profile experiment) using *.../Polynomials/s1_Poggensee/Main_polynomials_Pog1_mtp.m*. THese script generate the following .mat files

- muscle_spanning_joint_INFO.... (info on indexes of muscles that span each dof)
- MuscleInfo_... (coefficient of the polynomials)

### Save MT parameters

You also have to save the muscle-tendon paramaters of an opensim model to a .mat file. You can do this for example using the function *3dpredictsim\MuscleModel\Rajagopal2015\SaveMuscleParam.m*

### Create casadi functions

We then create casadi functions for the most important matlab functions used in the optimization. You can create these casadi functions using the script *CasadiFunctions_all_mtp_CreateDefault.m*. On Line 11 you can select the folder where the casadi functions will be saved.

```matlab
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp_Default';
```

Note that you if you want to make changes in the parameters of the hill model, you have to change this here. We might want to adjust this workflow if we run a sensitivity analysis. The advantage of this implementation is that you only have to create the casadifunctions once if you want to simulate different exoskeleton controllers (i.e. you can adapt the exoskeleton control in the next step).

### Run predictive simulation

You can use the function *../OCP/f_PredSim_PoggenSee2020.m* to run the predictive simulation. You have to provide a matlab structure as input argument for this function with the "settings" for this function. You can find a good example on how to use this function in  *../RunSim/Example_SimulateExo.m*. Note that this structure is updated in the function *f_PredSim_PoggenSee2020.m* with the default settings (i.e. settings not provided as input) on line 29. 

```matlab
S = GetDefaultSettings(S);
```

In this settings files, you have to point to the right folder with casadi function, .dll file for the skeleton dynamics, number of mesh points, desired movement speed, number of threads and so on.... (see details in example and in *GetDefaultSettings*).

 ### Post process results

The function *f_PredSim_PoggenSee2020* only saves the optimization variables (i.e. states and control). Using the function *f_LoadSim_PoggenSee2020_DefaultS* you can analyse your simulation. This analysis includes

- create a .mot file to visualise the simulation result in OpenSim
- computation metabolic energy & cost of transport
- compute spatio-temporal outputs
- ...

