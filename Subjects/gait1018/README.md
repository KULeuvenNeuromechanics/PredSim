Example 2D model with 10 degrees of freedom and 9 muscles per leg. 
Original model file can be found in `*/OpenSim 4.3/Resources/Models/Gait10dof18musc/gait10dof18musc.osim`.

Since the provided model is not compatible with PredSim, it was adapted using 
https://github.com/KULeuvenNeuromechanics/PredSim/blob/add_implementation_2D_model/AdaptOpenSimModel/AdaptOpenSimModel.m to include:
- torque actuator lumbar joint
- foot-ground contact spheres
- polynomial description of knee axis translation in function of knee angle (i.e. Yamaguchi knee)
