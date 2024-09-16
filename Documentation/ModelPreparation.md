
# Preparing a musculoskeletal model for use with PredSim

## What information does PredSim use to describe the musculoskeletal dynamics?

#### OpenSim model
- bodies, joints, coordinates
- inverse dynamics (via OpenSimAD)
- muscles (names, geometry, Hill-type params)
- ligaments (names, geometry, params)
- coordinate actuators
- contact elements
- coordinate bounds

#### Settings
- bounds on coordinates and activations
- scale muscle params (extended)
- coordinate stiffness and damping
- coordinate limit torques
- ligament stress-strain model (see [template model](../ModelComponents/ligamentForceLength_template.m))

#### Hardcoded
- muscle contraction and activation dynamics - F. De Groote, A. L. Kinney, A. V. Rao, and B. J. Fregly, “Evaluation of Direct Collocation Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem,” Annals of biomedical engineering, vol. 44, no. 10, pp. 2922–2936, 2016.
- muscle activation and deactivation time constants - [0.015s and 0.060s resp.]
- [ratio of slow twitch fibres for each muscle](../PreProcessing/getSlowTwitchRatios.m) - used for metabolic energy model
- [specific tension of each muscle](../PreProcessing/getSpecificTensions.m) - used to calculate muscle mass for metabolic energy model
