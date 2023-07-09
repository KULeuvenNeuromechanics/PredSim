---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

*Before submitting a bug report, read the error messages in MATLAB. Some errors are there to inform you about inconsistent inputs for certain settings.* 

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Which script did you run?
2. What settings did you use?
3. What kind of OpenSim model did you use? Link to generic model, attach specific model, or describe model (bodies, joints, coordinates) if it is confidential.

**Log file**
The log files will tell us what the error is, and where it occurred. 
When running main.m, you can find a log file of each simulation in the folder `S.subject.save_folder`. The name of this file ends with `_log.txt`. If there is no log file, copy the MATLAB command window to a .txt file and attach that instead.

If you run main.m and encountered an error before seeing "[...PreProcessing done. Time elapsed ** s](https://github.com/KULeuvenNeuromechanics/PredSim/blob/81068eae09c1f75a514ba0885df579e84dced472/Tests/Falisse_et_al_2022_Results/Falisse_et_al_2022_v1_log.txt#L16)" in the MATLAB command window, we need more detailed logs. Set `S.OpenSimADOptions.verbose_mode = true;` and run again. This will not solve your problem, but the new log file will give us more information.

**Software versions**
 - Windows 10
 - MATLAB ...
 - OpenSim ...
 - CasADi ...
 - Microsoft Visual Studio ...
 - CMake ...

**Additional context**
Add any other context about the problem here.
