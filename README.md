# neuronal_dynamics

Code for ["Between-area communication through the lens of within-area neuronal dynamics"](https://www.biorxiv.org/content/10.1101/2022.04.11.487906v3), Gozel and Doiron, *bioRxiv*, 2023.

## simulate

This folder contains all the scripts needed to run the simulations of the three-layer network of recurrently connected spiking neurons. In particular:
* Sim_2layers_cluster_sameParams.m : Run this script on a high-performance computing cluster to simulate the full network for several trials and several parameter sets (including the standard parameter set). (NB: Follow the instructions at the top of the script for set-up.)
* Sim_1layer_cluster_sameParams.m : Run this script on a high-performance computing cluster to simulate an input--output network used to perform linearity analysis. (NB: Follow the instructions at the top of the script for set-up.)
