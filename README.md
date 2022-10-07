# Neuromuscular controller parameter estimation

This repository estimates control parameters of a neuromuscular controller to track inverse dynamic moments during steady-state walking and in response to pelvis push and pull perturbations. More details in [ref to publication].

### Approach

We use an inverse skeleton dynamics approach. This means that we compute the inverse dynamic joint moments from motion capture data an try to track the ID moments with a model consisting of hill-type muscles driven by reflexes.

### Installation instructions

- Download and install casadi in matlab (https://web.casadi.org/)
- Install the opensim matlab API (https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab)
- clone this repository and add the folder (and subfolder) to your matlab path (e.g. using the command addpath(genpath(...)))
- Adapt in each script the main path (*Set.MainPath      ='C:\Users\mat950\Documents\Software\Sim\GeyerAnkle_ParamEstimation')* and point to the folder where you cloned this repository

### Download Experimental data

You can download the experimental data of the pelvis push and pull perturbations [1] using

- Run matlab function *GetDataVlutters2018.m* (by typing GetDataVlutters2018.m in your matlab command prompt)

### Code structure

- Adapt in each script the main path (*Set.MainPath      ='C:\Users\mat950\Documents\Software\Sim\GeyerAnkle_ParamEstimation')* and point to the folder where you cloned this repository
- First run the script: *Scripts/CreateCasadiFunctions_Shooting.m* to casadif functions for the system dynamics and the integration scheme. This function makes it easier (and faster) to build the optimization problem
- Estimate parameters using: *Scripts/ParamEst_PelvisPerturb_Shooting.m*. This script formulates and solves the optimizaton problem
- *PredictNovelPerturbationMagnitude.m*: runs a forward simulation with the identified gains on perturabtions that were not included in the parameter estimation process.
- *PlotResultsParamID.m*: Plots results
- *ParamEst_PelvisPerturb_Shooting_MinimalExample.m*. This scripts solve the optimization problem on a simple dataset (with a limited number of gait cycles)

### Background information

I use a matlab structure (names "Set") as input in my functions to do the flow control and set important parmaeters in the optimimzation. With this structure you can change the default settings (defined in ***GetDefaultSettings.m***). (Please read the function *GetDefaultSettings.m*  for all possible settings)

### References

[1] Vlutters, M., van Asseldonk, E.H.F. & van der Kooij, H. Lower extremity joint-level responses to pelvis perturbation during human walking. *Sci Rep* **8**, 14621 (2018). https://doi.org/10.1038/s41598-018-32839-8 