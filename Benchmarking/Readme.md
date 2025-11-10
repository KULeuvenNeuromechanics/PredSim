## Benchmarking

info about why benchmarking and refer to paper


## Content datafiles

Data is stored as a json file. This json file contains non-dimensionalised data
in the following structure.

- ik: matrix with inverse kinematics data for one gait cycle (nfr x ndof)
- ik_std: matrix with standard deviation inverse kinematics data for one gait cycle (nfr x ndof)
- ik_header: header inverse kinematics
- id: matrix with inverse dynamics data for one gait cycle (nfr x ndof)
- id_std: matrix with standard deviation inverse dynamics data for one gait cycle (nfr x ndof)
- id: header inverse dynamics
- grf_r: matrix with ground reaction force on right leg (nfr x 3)
- grf_r_std: matrix with std ground reaction force on right leg
- grf_l: matrix with ground reaction force on left leg (nfr x 3)
- grf_l_std: matrix with std ground reaction force on left leg (nfr x 3)
- grf_header: header GRF (default is {'Fx','Fy','Fz'};
- stride_frequency: stride frequency
- Pmetab_mean: mean metabolic power
- subject_height: (average) height of participant(s) (in meter !)
- subject_mass: (average) mass of the participant(s) (in kg !)
- speed: walking speed in m/s
- slope: slope (e.g. slope = 0.08 for walking on an 8% incline)
- added_mass:  total added mass to subject (in kg)
- location_added_mass: name of body (without _l or _r)
- study: name of the study (e.g. browning2008)

Feel free to add other outputs (like exoskeleton, exoskeleton_controller) to this. Everything you forget 
to add will be treated as empty.

Tot nondim outputs based on:
- frequency:  sqrt(g/l)
- moments:    m*g*l
- forces:     m*g
- Pmetab:     m*g^1.5*sqrt(l)

## identifier for simulation

We use an unique identifier for each experimental condition and the connected predictive simulation. This makes it easier to connect experimental data to simulation. This identifier is a string. If you add experimental data to () it would be nice to also add this identifier to the list below.



**Koelewijn 2019:**

- 0.8 m/s, 0% slope: **koelewijn2019_0p8ms**

- 1.3m/s, 0% slope: **koelewijn2019_1p3ms**

- 0.8, -8%slope: **koelewijn2019_0p8ms_8decline**

- 1.3, -8%slope: **koelewijn2019_1p3ms_8decline**

- 0.8, 8%slope: **koelewijn2019_0p8ms_8incline**

- 1.3, 8%slope: **koelewijn2019_1p3ms_8incline**


**Browning2008**

- xkg femur: browning2008_femurxkg
- xkg foot: browning2008_footxkg
- xkg pelvis: browning2008_pelvisxkg
- xkg tibia: browning2008_tibiaxkg

**VanDerZee2022**

- x speed (in m/s): vanderzee2022_xms 
- example: vanderzee2022_0p7ms for walking at 0.7 m/s


**Gomenuka2014**

- walking on slope x with added mass y at z m/s: gomenuka2014_slope_xpct_ykmh_mass_zpctmass
- example for walking on slope 7 percent at 3 km/h and 25% added mass: gomenuka_slope_7pct_3kmh_25pctmass


**Gait Speeds**
- walking at x m/s: gait_speeds_xms
- example: for walking at 0.85 ms: gait_speeds_0p85ms





### Idea benchmarking workflow

The main idea is that you can test your simulation workflow on a set of gait conditions and compare your simulations to experiments. 

**Input:**

- Simulation model

- Settings used for you simulations 

**Output**

- Simulation results

- experimental values associated with motion, preferably scaled for the person (non-dimentionalised for paper and re-computed for antropometry of specific model used in simulation)
  
  - spatio-temporal
  
  - kinematics
  
  - kinetics
  
  - metabolic power

**benchmark conditions:**

- gait at various speeds [Abe2025, VanDerZee?]

- walking on a slope [Koelewijn?, McDonald?]

- walking with added mass [schertzer dataset ?]
