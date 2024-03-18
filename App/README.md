
App developed for Kinderuniversiteit Leuven - Gasthuisberg 2022 and updated for a research internship.
=================

# Installation

- For the base installation, please refer to the README in the main folder.

Specifics for the app:
- in Vitruvian_Man_NL.mlap through the code view: in the function setPaths add your computer and the respective paths
- in PredSim_wrapper_for_app.m: on line 134, make sure you have the correct compiler selected for your system

After changing the .mlap file, you can export the app to a .m file; that is easier to use. 

Run the exported .m file and see if it works.

# Simulating

Specifics using the GUI
- Groep: give your group name, this will be the subfolder where results will be stored
- Modelnaam: this name should be unique for each model. I.e. if changeing model length or model mass you should give a new name.
