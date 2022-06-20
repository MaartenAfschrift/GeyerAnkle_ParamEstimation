# GeyerAnkle Parameter estimation

Objective: Estimate reflex parameter of a hill-type muscle model by tracking experimental data of (1) unperturbed walking, (2) treadmill belts speed perturabtions during walking and (3) pelvis push perturbation during walking.

### Approach

We use an inverse skeleton dynamics appraoch. This means that we compute the inverse dynamic joint moments from motion capture data an try to track the ID moments with a model consisting of hill-type muscles driven by reflexes.

### Code structure

*To Do*



## Background information



### Code structure

#### Settings:

I use a matlab structure (names "Set") as input in my functions to do the flow control and set important parmaeters in the optimimzation. With this structure you can change the default settings (defined in ***GetDefaultSettings.m***).

Fields:

- **BoolPerturb:** dataset contains perturbations?
- **COMfb:** true if you want to use COMfb
- **PercGaitCycle:** true if you want to compute the deviation in COM kinematics as a percentage of the gait cycle
- **COM_std**: Boolean of you want to use two times the standard deviation in unpeturbed com kinematics as a bound around the average COM trajectories when compute the COM deviation for the unperturbed trajectories.
- .... [Please read the function *GetDefaultSettings.m*  for all possible settings]



### Model info

#### Evalating muscle-tendon length using polynomials

We fitted polynomial function to compute the muscle-tendon length and moment arms as a function of the joint kinematics. You can find this polynomial fitting on a "training dataset" here: Functions\PolynomialFitting\Main_Polynomials.m

#### Inverse kinematics and dynamics

Inverse kinematics and kinetics were computed in OpenSim.



