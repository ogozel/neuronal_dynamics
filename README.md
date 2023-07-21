# neuronal_dynamics

Code for ["Between-area communication through the lens of within-area neuronal dynamics"](https://www.biorxiv.org/content/10.1101/2022.04.11.487906v3), Gozel and Doiron, *bioRxiv*, 2023.


## simulate

This folder contains all the scripts needed to run the simulations of the three-layer network of recurrently connected spiking neurons. In particular:
* Sim_2layers_cluster_sameParams.m : Run this script on a HPC cluster to simulate the full network for several trials and several parameter sets (including the standard parameter set). (NB: Follow the instructions at the top of the script for set-up.)
* Sim_1layer_cluster_sameParams.m : Run this script on a HPC cluster to simulate an input--output network used to perform linearity analysis. (NB: Follow the instructions at the top of the script for set-up.)
* runslurm_simLinear.sh : example of a bash script to run the simulations for the linearity analysis on a HPC cluster using a slurm scheduler
* get_combinedData.m : Run this script on the HPC cluster to combine the spiking data from the different trials with identical network parameters


## analyze

This folder contains the scripts needed to do the analyses.

* compute_FAcomSub.m : Run this script on a high-performance computing cluster to compute the within-area dimensionality of the Sender and Receiver networks, as well as their between-area communication. This scripts reads the data from the files with the combined spiking data of all the trials of each parameter set. (NB: follow the instructions at the top of the file)
* plot_FAcomSub.m : Once 'compute_FAcomSub.m' has been run and the results have been saved, this script can be run to plot many different analyses of the relationship of within-area dimensionality, prediction performance of the communication subspace, and radius of the disc from which neurons are sampled.
