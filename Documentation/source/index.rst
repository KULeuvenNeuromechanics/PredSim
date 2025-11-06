

Predictive Simulations of Locomotion
====================================


This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of locomotion. The simulation framework is developed by Falisse *et al.*. 
The implementation in this repo is aimed at letting you run simulations with your customized musculoskeletal models.
If you want to reproduce specific published results, see replicate results.

If you use PredSim for research, please cite 

* Falisse A, Serrancoli G, Dembia C, Gillis J, Jonkers J, De Groote F. 2019 Rapid predictive simulations with complex musculoskeletal models suggest that diverse healthy and pathological human gaits can emerge from similar control strategies. Journal of the Royal Society Interface 16: 20190402. http://dx.doi.org/10.1098/rsif.2019.0402.

* Lars D’Hondt, Antoine Falisse, Dhruv Gupta, Bram Van Den Bosch, Tom J. W. Buurke, Míriam Febrer-Nafría, Ines Vandekerckhove, Maarten Afschrift, and Friedl De Groote, PredSim: A Framework for Rapid Predictive Simulations of Locomotion, 10th IEEE RAS/EMBS International Conference for Biomedical Robotics and Biomechatronics (BioRob), Heidelberg, Germany, 2024, pp. 1208-1213, http://dx.doi.org/10.1109/BioRob60516.2024.10719735.


This repository is a work in progress.

⚠️ **Troubleshooting Issues**:

   If you encounter an error or get an unexpected output, let us know. Does the answer to your question/remark probably require code changes?
   
   - Yes -> use the `issues <https://github.com/KULeuvenNeuromechanics/PredSim/issues>`__. Please look at the existing issues before creating a new one; you might not be the first person to have this problem.
   - No -> use the `discussions <https://github.com/KULeuvenNeuromechanics/PredSim/discussions>`__. 

🚧 **Suggesting Improvements**:

   Feel free to suggest improvements. Submit a pull request with the changes (`See Guidelines for contributing <https://github.com/KULeuvenNeuromechanics/PredSim/wiki/Guidelines-for-contributing-code>`__),
   or post it in the discussion section on `Feature requests & Suggestions <https://github.com/KULeuvenNeuromechanics/PredSim/discussions/categories/feature-requests-suggestions>`__.

📯❔ **Sharing Best Practices and General Questions**:

   To share best practices, publications, or more general questions, we encourage you to use the discussions section.


Table of Contents
-----------------

.. toctree::
   :maxdepth: 1
   :caption: Setup

   Installing PredSim <InstallingPredSim.rst>
   Getting started <GettingStarted.rst>
   Run PredSim on a cluster <SetupCluster>
   Creating a private fork <PrivateForkPredSim>

.. toctree::
   :maxdepth: 1
   :caption: Documentation

   Settings <SettingsOverview.md>
   Using custom models <ModelPreparation.md>
   Accessing results <ResultsStructure.md>
   Structure of the code <https://github.com/KULeuvenNeuromechanics/PredSim/wiki/Structure-of-the-code>
   Troubleshooting <Troubleshooting.rst>
   Release notes <https://github.com/KULeuvenNeuromechanics/PredSim/releases>

.. toctree::
   :maxdepth: 1
   :caption: Community

   User forum <https://github.com/KULeuvenNeuromechanics/PredSim/discussions>
   Contributing <https://github.com/KULeuvenNeuromechanics/PredSim/wiki/Guidelines-for-contributing-code>
   Publications <PUBLICATIONS>

.. toctree::
   :maxdepth: 1
   :caption: Tutorials and examples

   Simple examples <https://github.com/KULeuvenNeuromechanics/PredSim/tree/master/Examples>
   Workshop September 2024 <https://github.com/KULeuvenNeuromechanics/PredSim-workshop-bcn-2024>
   
