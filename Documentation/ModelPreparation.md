
# Preparing a musculoskeletal model for use with PredSim

- [What information does PredSim use to describe the musculoskeletal dynamics?](#what-information-does-predsim-use-to-describe-the-musculoskeletal-dynamics)
- [Limitations on the OpenSim model](#limitations-on-the-opensim-model)
- [Using a very different model](#using-a-very-different-model)


## What information does PredSim use to describe the musculoskeletal dynamics?

#### OpenSim model
- bodies, joints, coordinates
- inverse dynamics (via OpenSimAD)
- muscles: names, geometry, Hill-type parameters (FMo, lMo, alphao, vMmax, lTs)
- ligaments: names, geometry, parameters (resting length, pcsa force)
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
- muscle activation and deactivation time constants - [0.015s and 0.060s resp.](../PreProcessing/read_and_scale_MTparameters.m)
- [ratio of slow twitch fibres for each muscle](../PreProcessing/getSlowTwitchRatios.m) - used for metabolic energy model
- [specific tension of each muscle](../PreProcessing/getSpecificTensions.m) - used to calculate muscle mass for metabolic energy model


## Limitations on the OpenSim model

- Your model should not have locked joints. Locked joints would technically require having kinematic constraints, which are not supported. Replace them with weld joints instead.
- Constraints on coordinates will be ignored (eg, coupling constraints).
- Using SimmSplines to describe coordinates (e.g. Yamaguchi knee model) is not supported as the implementation in OpenSim is not really compatible with algorithmic differentiation. Change them to Polynomials instead. [_AdaptOpenSimModel.m_](../AdaptOpenSimModel/AdaptOpenSimModel.m) takes care of changing SimmSplines present in joint definitions to polynomials. GeometryPaths can contain SimmSplines.
- Your model needs to have contact elements that interact with the ground. Only *SmoothSphereHalfSpaceForce* contact forces are supported. You can use [_AdaptOpenSimModel.m_](../AdaptOpenSimModel/AdaptOpenSimModel.m) to add contact geometries and forces to your model. You can also scale the radius, stiffness and dissipation of the contact spheres.
- Your model can have any Hill-type muscle model, but it will be implemented as a [DeGroote-Fregly muscle](https://doi.org/10.1007/s10439-016-1591-9).
- Torque/force actuators of the class *ActivationCoordinateActuator* are supported. You can add actuators by running [_AdaptOpenSimModel.m_](../AdaptOpenSimModel/AdaptOpenSimModel.m). Actuators are not required.


## Using a very different model
- non-human
- with different joint definitions than the example models (e.g. knee flexion is positive)

Required changes w.r.t. default settings
- S.bounds.default_coordinate_bounds: create different table, or set to empty to use coordinate bounds from opensim model
- S.misc.default_msk_geom_bounds: create different table, or set to empty to use coordinate bounds from opensim model
- S.subject.default_coord_lim_torq_coeff: create different table
- S.bounds.distanceConstraints: prevent limbs from clipping through eachother/torso
- S.subject.base_joints_legs and S.subject.base_joints_arms: identify limbs

Changes in hardcoded parts might be required
- muscle activation and deactivation time constants - [0.015s and 0.060s resp.](../PreProcessing/read_and_scale_MTparameters.m)
- [ratio of slow twitch fibres for each muscle](../PreProcessing/getSlowTwitchRatios.m) - used for metabolic energy model
- [specific tension of each muscle](../PreProcessing/getSpecificTensions.m) - used to calculate muscle mass for metabolic energy model
