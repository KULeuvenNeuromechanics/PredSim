
Getting started
================

.. _dependencies:
Required software
-----------------
To run this code you need to have the following software on your machine:


**MATLAB**

- version: R2019b - R2025b
- toolboxes:
    - `Statistics and Machine Learning Toolbox <https://nl.mathworks.com/products/statistics.html>`__
    - `Signal Processing Toolbox <https://nl.mathworks.com/products/signal.html>`__
    - `Parallel Computing Toolbox <https://nl.mathworks.com/products/parallel-computing.html>`__ (optional)

.. _setup_opensim:
**OpenSim**

- `download <https://simtk.org/projects/opensim>`__
- version: 4.3 or later. Older versions do not work.
- Get the OpenSim API running on MATLAB. See `Setting up your Matlab Scripting Environment
 <https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab#ScriptingwithMatlab-MatlabSetupSettingupyourMatlabScriptingEnvironment>`_.
 Make sure to follow the steps exactly as described.

**CasADi**

- `download <https://web.casadi.org/get/>`__
- version: The code has been tested on CasADi 3.5.5 and later.

**Microsoft Visual Studio**

- `download <https://visualstudio.microsoft.com/>`__
- version: The code has been tested on MSVS Community 2015, 2017, 2019, and 2022.
- In Visual Studio Installer, `select to include Desktop development with C++ <./_static/fig_MSVS.png>`__

.. _setup_cmake:
**CMake**

- `download <https://cmake.org/download/>`__
- version: The code has been tested on CMake 3.22.0. 
- Add CMake to your system Path.

**git**

- `download <https://git-scm.com/download/win>`__
- version: The code has been tested on Git 2.40.0.windows.1. 
- Add Git to your system Path.



How to setup the code
---------------------
1. Fork this repository to your GitHub account. If you want the fork to be private, follow `these steps <./PrivateForkPredSim.md>`_ instead.
2. Clone the fork to your machine. If you have a computer with restricted permissions, make sure you have permission to run executables from the selected folder (For computers with KU Leuven BioMed Group policies, this is C:\GBW_MyPrograms\).
   *Do not download the code as zip.*
3. Add the location where you installed CasADi to the MATLAB search path (``addpath(genpath('C:/path/to/casadi'))``), or use ``S.solver.CasADi_path``.
4. Make sure the opensimAD submodule is installed. If PredSim/opensimAD/ is empty, open git command prompt, go to PredSim/, and run ``git submodule update --init``.

After performing these steps, run the main script. (Expected run time is 40 minutes, depending on hardware.) If you don't receive any errors, your results should be the same as `the reference result <https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Tests/ReferenceResults/Falisse_et_al_2022>`_. If that is the case, you have successfully installed and set up the code. You are ready to do your own simulations.

How to use the code
-------------------
The code is written such that as a user you only have to interact with `main.m <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/main.m>`_.
All user-defined `settings <./SettingsOverview.md>`_ are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings.

Unspecified settings use their default values. Defaults are described in the overview of `settings <./SettingsOverview.md>`_. A selection of impactful defaults is given below:

- ...
- ...


The `examples folder <https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Examples>`_ contains scripts illustrating how to use settings to:

- Predict a gait pattern resulting from a limited number of muscle synergies.
- Predict gait with assistance from an ankle exoskeleton.
- Run a sensitivity analysis.
