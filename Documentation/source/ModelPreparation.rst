
Preparing a musculoskeletal model for use with PredSim
======================================================



What information does PredSim use to describe the musculoskeletal dynamics?
---------------------------------------------------------------------------

**OpenSim model**

- bodies, joints, coordinates
- inverse dynamics (via OpenSimAD)
- muscles: names, geometry, Hill-type parameters (FMo, lMo, alphao, vMmax, lTs)
- ligaments: names, geometry, parameters (resting length, pcsa force)
- coordinate actuators
- contact elements
- coordinate bounds

**Settings**

- bounds on coordinates and activations
- scale muscle params (extended)
- coordinate stiffness and damping
- coordinate limit torques
- ligament stress-strain model (see `template model <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/ModelComponents/ligamentForceLength_template.m>`__)

**Hardcoded**

- muscle contraction and activation dynamics - :ref:`De Groote et al, 2016 <[De Groote 2016]>`
- muscle activation and deactivation time constants - `0.015s and 0.060s resp. <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/read_and_scale_MTparameters.m>`__
- `ratio of slow twitch fibres for each muscle <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/getSlowTwitchRatios.m>`__ - used for metabolic energy model
- `specific tension of each muscle <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/getSpecificTensions.m>`__ - used to calculate muscle mass for metabolic energy model


Limitations on the OpenSim model
--------------------------------

- Your model should not have locked joints. Locked joints would technically require having kinematic constraints, which are not supported. Replace them with weld joints instead.
- Constraints on coordinates will be ignored (eg, coupling constraints).
- Using SimmSplines to describe coordinates (e.g. Yamaguchi knee model) is not supported as the implementation in OpenSim is not really compatible with algorithmic differentiation. Change them to Polynomials instead. `AdaptOpenSimModel.m`_ takes care of changing SimmSplines present in joint definitions to polynomials. GeometryPaths can contain SimmSplines.
- Your model needs to have contact elements that interact with the ground. Only *SmoothSphereHalfSpaceForce* contact forces are supported. You can use `AdaptOpenSimModel.m`_ to add contact geometries and forces to your model. You can also scale the radius, stiffness and dissipation of the contact spheres.
- Your model can have any Hill-type muscle model, but it will be implemented as a `DeGroote-Fregly muscle <[De Groote 2016]>`.
- Torque/force actuators of the class *ActivationCoordinateActuator* are supported. You can add actuators by running `AdaptOpenSimModel.m`_. Actuators are not required.

.. _AdaptOpenSimModel.m: https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/AdaptOpenSimModel/AdaptOpenSimModel.m



Using a very different model
----------------------------

- non-human
- with different joint definitions than the example models (e.g. knee flexion is positive)

Required changes w.r.t. default settings

- S.bounds.default_coordinate_bounds: create different table, or set to empty to use coordinate bounds from opensim model
- S.misc.default_msk_geom_bounds: create different table, or set to empty to use coordinate bounds from opensim model
- S.subject.default_coord_lim_torq_coeff: create different table
- S.bounds.distanceConstraints: prevent limbs from clipping through eachother/torso
- S.subject.base_joints_legs and S.subject.base_joints_arms: identify limbs

Changes in hardcoded parts might be required

- muscle activation and deactivation time constants - `0.015s and 0.060s resp. <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/read_and_scale_MTparameters.m>`__
- `ratio of slow twitch fibres for each muscle <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/getSlowTwitchRatios.m>`__ - used for metabolic energy model
- `specific tension of each muscle <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/PreProcessing/getSpecificTensions.m>`__ - used to calculate muscle mass for metabolic energy model

