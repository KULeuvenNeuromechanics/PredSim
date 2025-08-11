


## R: results

The results are organised in a struct *R* that contains matrices for the different variables. 
The dimensions of the matrices depend on the model and settings that were used. 
In the overview below, the dimensions are expressed in function of:
- mesh: the number of meshpoints used to express the result.
If *S.misc.gaitmotion_type* is *FullGaitCycle*, the results use *S.solver.N_meshes* meshpoints.
If *S.misc.gaitmotion_type* is *HalfGaitCycle*, the results use *2\*S.solver.N_meshes* meshpoints.
- coordinates: the number of coordinates.
- muscles: the number of muscles.
- ligaments: the number of ligaments.
- actuators: the number of (coordinate) torque actuators.

The order of *mesh* represents 1 gait cycle, starting at the right side initial contact.
The order of *coordinates*, *muscles*, and *ligaments* is given in *R.colheaders*.


### R.S: settings used to generate these results
See the [overview of settings](../README.md#Required-Settings) for more information about the settings.

### R.objective
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| absoluteValues | values of the objective function terms, evaluated at the solution | 1 x 10 | [a.u.] |
| relativeValues | relative contribution of the objective function terms | 1 x 10 | [%] |
| labels | names of the objective function terms | 1 x 8 | [text] |

### R.time
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| mesh | timestamps of mesh points of simulation | 1 x (*S.solver.N_meshes* + 1) | [s] |
| coll | timestamps of mesh and collocation points of simulation | 1 x (*S.solver.N_meshes*\*(collocation_order + 1) + 1) | [s] |
| mesh_GC | timestamps of mesh points of gait cycle | 1 x (mesh + 1) | [s] |

### R.colheaders
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| coordinates | names of the coordinates, in the order that they are stored in the matrices with results | 1 x coordinates | [text] |
| muscles | names of the muscles, in the order that they are stored in the matrices with results | 1 x muscles | [text] |
| objective | names of the objective terms, in the order that they are stored in the matrices with results | 1 x 10 | [text] |
| GRF | names of the ground reaction forces directions, in the order that they are stored in the matrices with results | 1 x 3 | [text] |
| ligaments | names of the ligaments, in the order that they are stored in the matrices with results | 1 x ligaments | [text] |

### R.kinematics
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| Qs | coordinate positions | mesh x coordinates | [° or m] |
| Qdots | coordinate velocities | mesh x coordinates | [°/s or m/s] |
| Qdotdots | coordinate accelerations | mesh x coordinates | [°/s^2 or m/s^2] |
| Qs_rad | coordinate positions | mesh x coordinates | [rad or m] |
| Qdots_rad | coordinate velocities | mesh x coordinates | [rad/s or m/s] |
| Qdotdots_rad | coordinate accelerations | mesh x coordinates | [rad/s^2 or m/s^2] |

### R.muscles
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| a | muscle activation | mesh x muscles | [-] |
| da | time derivative of activation | mesh x muscles | [1/s] |
| FTtilde | force at the tendon (normalised to FMo) | mesh x muscles | [-] |
| dFTtilde | time derivative of force at the tendon (normalised to FMo) | mesh x muscles | [1/s] |
| e | muscle excitation | mesh x muscles | [-] |
| lMT | muscle-tendon length | mesh x muscles | [m] |
| vMT | muscle-tendon lengthening velocity | mesh x muscles | [m/s] |
| dM | muscle-tendon moment arm | mesh x muscles x coordinates | [m] |
| FT | force at the tendon  | mesh x muscles | [N] |
| Fce | force of the contractile element (CE)  | mesh x muscles | [N] |
| Fpass | force of the parallel elasticity (PE)  | mesh x muscles | [N] |
| Fiso | max isometric force of the contractile element (CE) (normalised to FMo) | mesh x muscles | [-] |
| lM | muscle fibre length | mesh x muscles | [m] |
| lMtilde | muscle fibre length (normalised to lMo) | mesh x muscles | [-] |
| vM | muscle fibre lengthening velocity | mesh x muscles | [m/s] |
| vMtilde | muscle fibre lengthening velocity (normalised to vMmax) | mesh x muscles | [1/s] |
| lT | tendon length | mesh x muscles | [m] |
| vT | tendon lengthening velocity | mesh x muscles | [m/s] |

Optional, if muscle synergies are implemented:
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| SynH_r | synergy activation, right muscles | mesh x synergies_right | [-] |
| SynH_l | synergy activation, left muscles | mesh x synergies_left | [-] |
| SynW_r | synergy weights, right muscles | synergies_right x muscles_right | [-] |
| SynW_l | synergy weights, left muscles | synergies_left x muscles_left | [-] |

<img src="./Hill-type%20muscle%20model.svg" width="400" height="auto">

Hill-type muscle



### R.torque_actuators
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| a | actuator activation | mesh x actuators | [-] |
| e | actuator excitation | mesh x actuators | [-] |
| T | actuator torque | mesh x actuators | [Nm] |


### R.ground_reaction
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| threshold | threshold for the vertical ground reaction force to detect stance | 1 x 1 | [N] |
| initial_contact_side | side where initial contact was detected before post-processing | 1 x 1 | [text] |
| GRF_r | ground reaction forces right side | mesh x 3 | [N] |
| GRF_l | ground reaction forces left side | mesh x 3 | [N] |
| GRM_r | ground reaction moments right side | mesh x 3 | [Nm] |
| GRM_l | ground reaction moments left side | mesh x 3 | [Nm] |
| COP_r | centre of pressure right side (in ground frame) | mesh x 3 | [m] |
| COP_l | centre of pressure left side (in ground frame) | mesh x 3 | [m] |
| idx_stance_r | indices of meshpoints with ground contact on right side | <mesh x 1 | [-] |
| idx_stance_l | indices of meshpoints with ground contact on left side | <mesh x 1 | [-] |
| idx_double_support | indices of meshpoints with ground contact on both sides | <mesh x 1 | [-] |

### R.misc
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| body_mass | body mass | 1 x 1 | [kg] |
| body_weight | body weight | 1 x 1 | [N] |

### R.kinetics
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| T_ID | inverse dynamics (inertial and contact) | mesh x coordinates | [Nm or N] |
| T_spring | coordinate moments by springs | mesh x coordinates | [Nm or N] |
| T_damping | coordinate moments by dampers | mesh x coordinates | [Nm or N] |
| T_limit | coordinate moments by limit torques | mesh x coordinates | [Nm or N] |

### R.ligaments
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| l | ligament length | mesh x ligaments | [m] |
| v | ligament lengthening velocity | mesh x ligaments | [m/s] |
| F | ligament force | mesh x ligaments | [N] |
| strain | ligament strain | mesh x ligaments | [%] |
| stress | ligament stress | mesh x ligaments | [MPa] |
| moment | coordinate moment (or force) due to ligaments | mesh x coordinates | [Nm or N] |
| power | coordinate power due to ligaments | mesh x coordinates | [W] |

### R.spatiotemp
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| stance_r | duration of stance phase of right leg, relative to gait cycle | 1 x 1 | [%] |
| swing_r | duration of swing phase of right leg, relative to gait cycle | 1 x 1 | [%] |
| stance_l | duration of stance phase of left leg, relative to gait cycle | 1 x 1 | [%] |
| swing_l | duration of swing phase of left leg, relative to gait cycle | 1 x 1 | [%] |
| double_support | duration of double support phase, relative to gait cycle | 1 x 1 | [%] |
| stride_freq | stride frequency | 1 x 1 | [Hz] |
| step_width_COP | step width, calculated based on centre of pressure | 1 x 1 | [m] |
| dist_trav | distance traveled, calculated based on floating base forward position | 1 x 1 | [m] |

### R.metabolics

#### R.metabolics.Bhargava2004

| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| Edot_gait | Metabolic energy rate of each muscle. $\dot{E}$ [in equation 1](https://doi.org/10.1016/S0021-9290(03)00239-2) | mesh x muscles | [W/kg] |
| Adot | Activation heat rate of each muscle. $\dot{A}$ [in equation 3](https://doi.org/10.1016/S0021-9290(03)00239-2) | mesh x muscles | [W/kg] |
| Mdot | Maintenance heat rate of each muscle. $\dot{M}$ [in equation 7](https://doi.org/10.1016/S0021-9290(03)00239-2) | mesh x muscles | [W/kg] |
| Sdot | Shortening heat rate of each muscle. $\dot{S}$ [in equation 8](https://doi.org/10.1016/S0021-9290(03)00239-2) | mesh x muscles | [W/kg] |
| Wdot | Work rate of each muscle. $\dot{W}$ [in equation 12](https://doi.org/10.1016/S0021-9290(03)00239-2) | mesh x muscles | [W/kg] |
| Edot_incl_basal | Total metabolic energy rate. | mesh x 1 | [W/kg] |
| COT | Cost of transport  | 1 x 1 | [J/kg/m] |


