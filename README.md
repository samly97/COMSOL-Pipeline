# COMSOL-Pipeline

Automatically generated 2D microstructures and study of Lithium-ion batteries.

# Instructions

1. Start up a COMSOL Multiphysics Server instance.
2. Run `comsol_init.m` to connect. Note that the COMSOL_DIR will be different depending where COMSOL is installed on your computer.
3. Change the `PATH` and `RESULTS_PATH` variables in `main.m` corresponding to where this code lives and where you want to save the results.
4. The minimum and maximum porosities could be modified in lines `6-7`. The number of microstructures can be modified on line `10`. Specify more C-rates in line `34`. Other settings could be modified to your desires.

It is recommend to run one simulation to see what materials are modeled and under which settings. Users could model other materials (electrolytes, active materials, etc.), but the code in `comsol_fns.m` would have to be modified. In that case, users could use export their COMSOL simulations to MATLAB code and refactor as required.
