
Getting started
===============


The code is written such that as a user you only have to interact with `main.m <https://github.com/KULeuvenNeuromechanics/PredSim/blob/master/main.m>`__, or a similar script.
All user-defined `settings <./SettingsOverview.md>`__ are stored in structure *S*. In main.m you have to specify the required settings and are free to change/add the optional settings.

**1. Initialise settings**

.. code-block::
    :caption: Initialise base settings

    [S] = initializeSettings();


.. code-block::
    :caption: Initialise base settings and model-specific defaults

    [S] = initializeSettings('Falisse_et_al_2022');


**2. Add more settings**

.. code-block::

    S.subject.name = 'Falisse_et_al_2022';

.. code-block::

    S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name);

.. code-block::

    S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
    S.solver.IG_selection_gaitCyclePercent = 100;

.. code-block::

    S.solver.IG_selection = 'quasi-random';



A selection of impactful defaults is given below:

- ...
- ...

.. code-block::

    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

**3. Start the simulation**

.. code-block::

    runPredSim(S, osim_path)

Unspecified settings use their default values. Defaults are described in the overview of `settings <./SettingsOverview.md>`__. 




**More examples**

The `examples folder <https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Examples>`__ contains scripts illustrating how to use settings to:

- Predict a gait pattern resulting from a limited number of muscle synergies.
- Predict gait with assistance from an ankle exoskeleton.
- Run a sensitivity analysis.
