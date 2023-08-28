Vitruvian Man
=================

- Original app developed for Kinderuniversiteit Leuven - Gasthuisberg 2022 (Vitruvian_Man_NL_exported).
- App translated and updated for TGCS 2023 (Vitruvian_Man_EN_exported_ratios).

## Installation

1. For the base installation, please refer to [the README in the main folder](../README.md).
2. Update [getPaths.m](./getPaths.m) with the relevant paths for your computer.
3. (Optional) For realtime highscores across multiple computers, you can write results to a google form via [autoFill.py](./autoFill.py). This needs additional setup:
    - have Python installed and working with MATLAB. See https://www.mathworks.com/products/matlab/matlab-and-python.html
    - for Python, install the library requests. This can be done through the command prompt and the command pip install requests
4. Run the "exported" file and see if it works.

If you get an error in the matlab command window and the last output before the error is "Creating new external function...", manually selecting the compiler might solve this. 
In [PredSim_wrapper_for_app](./PredSim_wrapper_for_app.m), uncomment the line with `S.OpenSimADOptions.compiler = ` and pass the compiler string for your visual studio installation. See [the README in the main folder](../README.md) for more info.



