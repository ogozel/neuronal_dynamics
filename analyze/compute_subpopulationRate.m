% Get subpopulation rate (spikes/s) as a function of time

% Run this script on the HPC cluster. Specify the data path accrodingly and
% load the file of interest. Then copy the saved file on local machine
% before running 'plot_subpopulationRateDifference_perDisc.m'


clear all
close all
clc


% Data path on the HPC cluster
dataFolder = '/user_data/ogozel/';


% seed1000_V2_sigma0d1 = load([dataFolder,'RecFeed2D_L1drivesL2_Uncorr_fixW1_L1_L2_sameParams_L2_sigmaRRFromI0d1_21s_seed1000_V2.mat']);
seed1000_V1_sigma0d3 = load([dataFolder,'RecFeed2D_L1drivesL2_Uncorr_fixW1_L1_L2_sameParams_L2_sigmaRRFromI0d3_21s_seed1000_V1.mat']);
% seed1000_V2_tau24 = load([dataFolder,'RecFeed2D_L1drivesL2_Uncorr_fixW1_L1_L2_sameParams_L2_tauIdecay24_21s_seed1000_V2.mat']);


%%% Choose the proper filename! %%%
filename = 'seed1000_V1_sigma0d3';
dataSpk = seed1000_V1_sigma0d3.s2;

t_start = 1000; % ms
t_end = 21000; % ms
bool_plot = 0;

% Total number of excitatory neurons
NE = 40000;

% Take all excitatory neurons of the grid
idxSelNeur = 1:NE;
    
% Compute subpopulation activity
[popSpk, timeBins] = fct_subpopulationActivity(dataSpk, t_start, t_end, idxSelNeur, bool_plot);

% Save results
filenameSave = [dataFolder, 'popRate_', filename, '_s2_1to21s'];
save(filenameSave, 'timeBins', 't_start', 't_end', 'popSpk', 'idxSelNeur', '-v7.3')








