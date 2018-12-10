# asynchronous_discontinuous_collective_behaviour
This repository contains Fortran 95 code for the model described in the paper "Disentangling and modeling interactions in fish with burst-and-coast swimming reveal distinct alignment and attraction behaviors" by Calovi et al. 2018. This is NOT the code released by the authors of the paper. Instead, it is a part my summer project with them which I continued into my Master's Thesis

# Description of the subfolders in the repository
1. Trajectory : Contains the f95 scripts required for obtaining the simulated trajectory data. This script can also be used to generate other data as mentioned in the description in the script.
2. Statistics : Contains two main f95 scripts. One is to evaluate the transient time of the dynamics of the sytem. The other, is to evaluate the group polarisation. Along with this,fish_func.f95(similar to the one in the Trajectory subfolder), which contains the necessary functions for the main scripts to compile and run, is also included. The remainining scripts in this subfolder can be used to clean the simulated data and plot them. 

# Notes
1. The file names are titled with a personal bias to be in sync with the namin in personal workspaces and other repositories. Do not make any correlation between the file names and what the function does unless explicitly mentioned.
