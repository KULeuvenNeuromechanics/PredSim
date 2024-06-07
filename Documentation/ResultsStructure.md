


## R | results

### R.S | settings used to generate these results

### R.objective
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
|absoluteValues| values of the objective function terms, evaluated at the solution | 1 x 8 | [a.u.] |
| relativeValues | relative contribution of the objective function terms | 1 x 8 | [%] |
| labels | names of the objective function terms | 1 x 8 | [text] |

### R.time
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| mesh | timestamps of mesh points of simulation | 1 x (mesh + 1) | [s] |
| coll | timestamps of mesh and collocation points of simulation | 1 x (mesh*collocation_order + 1) | [s] |
| mesh_GC | timestamps of mesh points of gait cycle | 1 x mesh[^1] | [s] |

### R.colheaders
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| coordinates | names of the coordinates, in the order that they are stored in the matrices with results | 1 x coordinates | [text] |
| muscles | names of the muscles, in the order that they are stored in the matrices with results | 1 x muscles | [text] |
| objective | names of the objective terms, in the order that they are stored in the matrices with results | 1 x 8 | [text] |
| GRF | names of the ground reaction forces directions, in the order that they are stored in the matrices with results | 1 x 3 | [text] |

### R.kinematics
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| Qs | coordinate positions | mesh[^1] x coordinates | [° or m] |
| Qdots | coordinate velocities | mesh[^1] x coordinates | [°/s or m/s] |
| Qdotdots | coordinate accelerations | mesh[^1] x coordinates | [°/s^2 or m/s^2] |
| Qs_rad | coordinate positions | mesh[^1] x coordinates | [rad or m] |
| Qdots_rad | coordinate velocities | mesh[^1] x coordinates | [rad/s or m/s] |
| Qdotdots_rad | coordinate accelerations | mesh[^1] x coordinates | [rad/s^2 or m/s^2] |

### R.muscles
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| a | muscle activation | mesh[^1] x muscles | [-] |
| da | time derivative of activation | mesh[^1] x muscles | [1/s] |
| FTtilde | force at the tendon (normalised to FMo) | mesh[^1] x muscles | [-] |
| dFTtilde | time derivative of force at the tendon (normalised to FMo) | mesh[^1] x muscles | [1/s] |
| e | muscle excitation | mesh[^1] x muscles | [-] |
| lMT | muscle-tendon length | mesh[^1] x muscles | [m] |
| vMT | muscle-tendon lengthening velocity | mesh[^1] x muscles | [m/s] |
| dM | muscle-tendon length | mesh[^1] x muscles x coordinates | [m] |
| FT | force at the tendon  | mesh[^1] x muscles | [N] |
| Fce | force of the contractile element (CE)  | mesh[^1] x muscles | [N] |
| Fpass | force of the parallel elasticity (PE)  | mesh[^1] x muscles | [N] |
| Fiso | max isometric force of the contractile element (CE) (normalised to FMo) | mesh[^1] x muscles | [-] |
| lM | muscle fibre length | mesh[^1] x muscles | [m] |
| lMtilde | muscle fibre length (normalised to lMo) | mesh[^1] x muscles | [-] |
| vM | muscle fibre lengthening velocity | mesh[^1] x muscles | [m/s] |
| vMtilde | muscle fibre lengthening velocity (normalised to vMmax) | mesh[^1] x muscles | [1/s] |
| lT | tendon length | mesh[^1] x muscles | [m] |
| vT | tendon lengthening velocity | mesh[^1] x muscles | [m/s] |



<img src="./Hill-type%20muscle%20model.svg" width="400" height="auto">

Hill-type muscle



### R.torque_actuators
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| a | actuator activation | mesh[^1] x actuators | [-] |
| e | actuator excitation | mesh[^1] x actuators | [-] |
| T | actuator torque | mesh[^1] x actuators | [-] |


### R.ground_reaction
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| GRF_r | ground reaction forces right side | mesh[^1] x 3 | [N] |
| GRF_l | ground reaction forces left side | mesh[^1] x 3 | [N] |
| GRM_r | ground reaction moments right side | mesh[^1] x 3 | [Nm] |
| GRM_l | ground reaction moments left side | mesh[^1] x 3 | [Nm] |
| COP_r | centre of pressure right side (in ground frame) | mesh[^1] x 3 | [m] |
| COP_l | centre of pressure left side (in ground frame) | mesh[^1] x 3 | [m] |
| idx_stance_r | indices of meshpoints with ground contact on right side | <mesh[^1] x 1 | [-] |
| idx_stance_l | indices of meshpoints with ground contact on left side | <mesh[^1] x 1 | [-] |

### R.misc

### R.kinetics
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| T_ID | inverse dynamics (inertial and contact) | mesh[^1] x coordinates | [Nm or N] |
| T_spring | coordinate moments by springs | mesh[^1] x coordinates | [Nm or N] |

### R.spatiotemp

### R.metabolics

[^1]: 