% Combining the data from 30 simulation trials
% This has to be run on a high-performance computing cluster.
% NB: modify the path to the folder, name of the files to read from and
% name of the file to save according to needs (lines 16-23)

% Load data where the spike counts have been saved (option.saveSpkCounts=1)
% E1 are excitatory spike counts in L1, I1 are inhibitory spike counts in L1
% E2 are excitatory spike counts in L2, I2 are inhibitory spike counts in L2
% All have format N x p:
%       N : nbr of neurons in the corresponding population
%       p : nbr of 50ms timewindows in the simulations

%% Modify according to needs

% Data path
dataPath = '../../user_data/ogozel/';

% Filename to read from
f_read = 'RecFeed2D_Uncorr_saveSpkCnts1_fixW1_L1_L2_sameParams_L2_sigmaRRFromI0d3';
% RecFeed2D_Uncorr_saveSpkCnts1_fixW1_L1_L2_sameParams_L1_tauIdecay24

% Name of the saved file with combined data
f_save = '30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L2_sigmaI0d3.mat';
% 30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L1_tauIdecay24.mat'


%% Combine the data from different simulation trials

% Intialize the matrices
Matrix_exc1 = [];
Matrix_inh1 = [];
Matrix_exc2 = [];
Matrix_inh2 = [];
Matrix_C = [];
Matrix_Cbar = [];
Matrix_COVd = [];
Matrix_COVbar = [];

nTrials = 30;
Tburn = 1; % 1 [s]
timewindow = 0.05; % 50 [ms]

% Loop over all trials with the same parameter set
for i=1:nTrials
    
    load([dataPath, f_read, '_ID', num2str(i), '.mat'])
    
    Matrix_exc1 = cat(2, Matrix_exc1, E1(:,(Tburn/timewindow)+1:end));
    Matrix_inh1 = cat(2, Matrix_inh1, I1(:,(Tburn/timewindow)+1:end));
    Matrix_exc2 = cat(2, Matrix_exc2, E2(:,(Tburn/timewindow)+1:end));
    Matrix_inh2 = cat(2, Matrix_inh2, I2(:,(Tburn/timewindow)+1:end));
    
    Matrix_C = cat(3, Matrix_C, C);
    Matrix_Cbar = cat(1, Matrix_Cbar, Cbar);
    Matrix_COVd = cat(3, Matrix_COVd, COV_d);
    Matrix_COVbar = cat(1, Matrix_COVbar, COVbar);
end

% Save the results on the high-performance computing cluster
save(dataPath, f_save, 'Matrix_exc1','Matrix_exc2','Matrix_inh1','Matrix_inh2','Matrix_C','Matrix_Cbar','Matrix_COVd','Matrix_COVbar')
