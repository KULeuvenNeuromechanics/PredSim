
Troubleshooting
===============

Check the `user forum`_ and `known bugs`_.

.. _user forum: https://github.com/KULeuvenNeuromechanics/PredSim/discussions
.. _known bugs: https://github.com/KULeuvenNeuromechanics/PredSim/issues?q=is%3Aissue%20state%3Aopen%20label%3Abug


Common problems and how to solve them
-------------------------------------

.. _err no opensimad:

.. dropdown:: The folder ``PredSim/opensimAD/`` is empty.

    1. Open git command prompt
    2. go to ``PredSim/``
    3. run ``git submodule update --init``.


.. _err no cmake:

.. dropdown:: A matlab error preceded by 
    ``'cmake' is not recognized as an internal or external command, operable program or batch file.``

    CMake is not :ref:`installed correctly  <setup_cmake>`.







