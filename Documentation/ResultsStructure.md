


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
| coll | timestamps of mesh and collocation points of simulation | 1 x (mesh*collocatio_order + 1) | [s] |
| mesh_GC | timestamps of mesh points of gait cycle | 1 x mesh [^1] | [s] |

### R.colheaders
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| coordinates | names of the coordinates, in the order that they are stored in the matrices with results | [text] |
| muscles | names of the muscles, in the order that they are stored in the matrices with results | [text] |
| objective | names of the objective terms, in the order that they are stored in the matrices with results | [text] |
| GRF | names of the ground reaction forces directions, in the order that they are stored in the matrices with results | [text] |

### R.kinematics
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| Qs | coordinate positions | [° or m] |
| Qdots | coordinate velocities | [°/s or m/s] |
| Qdotdots | coordinate accelerations | [°/s^2 or m/s^2] |
| Qs_rad | coordinate positions | [rad or m] |
| Qdots_rad | coordinate velocities | [rad/s or m/s] |
| Qdotdots_rad | coordinate accelerations | [rad/s^2 or m/s^2] |

### R.muscles
| Field | Description | Dimension | Unit |
|------ | ----------- | ---- | ---- |
| a | activation | [-] |
| da | time derivative of activation | [1/s] |
| FTtilde | force at the tendon (normalised to FMo) | [-] |
| dFTtilde | time derivative of force at the tendon (normalised to FMo) | [1/s] |
| e | excitation | [-] |
| lMT | muscle-tendon length | [m] |
| vMT | muscle-tendon lengthening velocity | [m/s] |
| dM | muscle-tendon length | [m] |


<img src="./Hill-type%20muscle%20model.svg" width="400" height="auto">

Hill-type muscle



### R.torque_actuators

### R.ground_reaction

### R.misc

### R.kinetics

### R.spatiotemp

### R.metabolics

[^1]: 