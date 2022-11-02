% Compute communication subspace between Sender and Receiver by randomly
% selecting neurons within discs of a specified radius (selected neurons
% should have mean firing rate > 2 Hz) in each area, as well as the
% within-area shared participation ratio.

% Run this script on cluster, and save results locally.
% Then plot results by running 'plot_FAcomSub.m' on the local machine

clear
close all
clc

bool_cluster = 1;

if bool_cluster
    % To perform FA to get the within-area shared covariance matrix
    addpath('../fa_Yu');
    % To compute the communication subspace
    addpath(genpath('../comSub_Semedo'));
    
    rng('shuffle');
    AI = getenv('SLURM_ARRAY_TASK_ID');
    data_folder='/user_data/ogozel/'; % On cluster
    
    bool_plot = 0; % no plotting on the cluster
else
    % To perform FA to get the within-area shared covariance matrix
    addpath('C:\Users\olivi\Dropbox\Projects\fa_Yu');
    % To compute the communication subspace
    addpath(genpath('C:\Users\olivi\Dropbox\Projects\communication-subspace'));
    
    data_folder = '../data_analysis/'; % On local machine
    
    bool_plot = 1; % plotting when on local machine!
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binDuration = 0.05; % each bin is 50ms
NE = 40000; % total number of excitatory neurons


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to choose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose which parameter we are interested in
thisParam = 'sigmaRRfromIinL2';
% 'sigmaRRfromIinL1'; 'sigmaRRfromIinL2';
% 'tauIdecayinL1'; 'tauIdecayinL2';

bool_print = 0; % print output on the command line

% Radius of the population of neurons we select
% NB: disc with radius=5 contains 81 neurons
Master_pop_radius = [5, 10:10:40, 100];

% Number of randomly selected neurons within the neuronal disc
num_selNeur = 50;

% Number of draws of randomly selected neurons for a given population radius
numDraws = 5;

% Force a low-dimensional representation
bool_forceLow = 0; % if 0: do not force a lowD representation
numDimForced = 5;

% Parameters for FA
zDimList = 1:num_selNeur;
numFolds = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(thisParam,'sigmaRRfromIinL1')
    
    TypePrefix = '30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L1_sigmaI';
    theseParams = [0.1:0.025:0.3];
    
    paramName = '\sigma_{I}';
    paramUnit = '';
    colororder = winter(length(theseParams));
    colormapName = 'winter';
    
elseif strcmp(thisParam,'sigmaRRfromIinL2')
    
    TypePrefix = '30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L2_sigmaI';
    theseParams = [0.1:0.025:0.3];
    
    paramName = '\sigma_{I}';
    paramUnit = '';
    colororder = winter(length(theseParams));
    colormapName = 'winter';
    
elseif strcmp(thisParam,'tauIdecayinL1')
    
    TypePrefix = '30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L1_tauIdecay';
    theseParams = [8:2:24];
    
    paramName = '\tau_{Id}';
    paramUnit = 'ms';
    colororder = copper(length(theseParams));
    colormapName = 'copper';
    
elseif strcmp(thisParam,'tauIdecayinL2')
    
    TypePrefix = '30trials_spikecountsAfterTburn_bin50ms_PoissonInput_sameParams_L2_tauIdecay';
    theseParams = [8:2:24];
    
    paramName = '\tau_{Id}';
    paramUnit = 'ms';
    colororder = copper(length(theseParams));
    colormapName = 'copper';
    
end



for p=1:length(theseParams)
    
    Master_idxSelected = cell(length(Master_pop_radius),numDraws,2); % indices of the selected neurons
    Master_meanFR = cell(length(Master_pop_radius),numDraws,2); % mean firing rate of the selected neurons
    
    % shared dimensionality
    Master_sharedPR = NaN(length(Master_pop_radius),numDraws,2); % participation ratio
    Master_shared95 = NaN(length(Master_pop_radius),numDraws,2); % number of dims to explain 95% of the shared variance
    
    Master_Cshared = cell(length(Master_pop_radius),numDraws,2); % shared covariance matrix
    Master_W = cell(length(Master_pop_radius),numDraws,2); % Cshared = WW^T
    Master_tildeU = cell(length(Master_pop_radius),numDraws,2); % eigenvectors of shared covariance matrix
    Master_tildeS2 = cell(length(Master_pop_radius),numDraws,2); % eigenvalues of shared covariance matrix
    Master_Cprivate  = cell(length(Master_pop_radius),numDraws,2); % private covariance matrix (diagonal elements)
    
    Master_ALLcomSub = NaN(length(Master_pop_radius),numDraws,4); % all communication subspaces, for several dimensionalities M
    Master_FinalComSub = cell(length(Master_pop_radius),numDraws,4); % communication subspace with the optimal dimensionality
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    filename = strrep(sprintf([TypePrefix,'%.03g'],theseParams(p)),'.','d')
    fname = [data_folder,filename,'.mat'];
    DataStruct = load(fname);
    
    % Size of selected populations
    pop_sizes = NaN(1,length(Master_pop_radius));
    
    for s=1:length(Master_pop_radius)
        
        % Current radius of the population of neurons we select
        pop_radius = Master_pop_radius(s);
        disp(['Population radius is ',num2str(pop_radius)])
        
        for d=1:numDraws
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Randomly select neurons
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Select a random disc of neurons (taking care of the periodic boundary conditions)
            center_ID = randperm(NE,1);
            idxDisc = fct_select_discNeurons(center_ID,pop_radius,sqrt(NE),sqrt(NE),bool_plot);
            
            % Size of the disc of (pre)selected neurons
            pop_sizes(s) = length(idxDisc);
            
            for thisLayer=1:2
                
                if thisLayer==1
                    fulldata = double(DataStruct.Matrix_exc1);
                elseif thisLayer==2
                    fulldata = double(DataStruct.Matrix_exc2);
                else
                    error('Problem: there are only 2 layers!')
                end
                
                % Transform the data from number of spikes/bin to number of spikes/s to get firing rate in [Hz]
                firingRate = fulldata/binDuration;
                meanFRperNeur = mean(firingRate,2); % average firing rate of each neuron [Hz]
                idxFRbigger = find(meanFRperNeur>2); % all neurons whose mean firing rate is >2Hz
                
                % Neurons of the selected disc whose mean firing rate is >2Hz
                idxDiscFRbigger = intersect(idxDisc,idxFRbigger);
                if length(idxDiscFRbigger) >= num_selNeur % Randomly select num_selNeur neurons from the subset
                    idxToKeep_wrtSubset = randperm(length(idxDiscFRbigger),num_selNeur);
                    idxToKeep = idxDiscFRbigger(idxToKeep_wrtSubset);
                else % in case there is less than 'num_selNeur' neurons that have the properties we want, keep all that do
                    idxToKeep = idxDiscFRbigger;
                end
                data = fulldata(idxToKeep,:);
                
                if thisLayer==1
                    data_sender = data;
                elseif thisLayer==2
                    data_receiver = data;
                end
                
                % Perform FA and compute within-area SHARED dimensionality (based on shared
                % covarariance matrix)
                [d_shared_PR, d_shared_95, eigenvectors, eigenvalues, sharedCovarianceMatrix, privateCovarianceElements, W] = ...
                    fct_compute_sharedDim(data,zDimList,numFolds,bool_forceLow,numDimForced,bool_print);
                
                % Save data
                Master_idxSelected(s,d,thisLayer) = {idxToKeep};
                Master_meanFR(s,d,thisLayer) = {meanFRperNeur(idxToKeep)};
                Master_sharedPR(s,d,thisLayer) = d_shared_PR;
                Master_shared95(s,d,thisLayer) = d_shared_95;
                Master_Cshared(s,d,thisLayer) = {sharedCovarianceMatrix};
                Master_W(s,d,thisLayer) = {W};
                Master_tildeU(s,d,thisLayer) = {eigenvectors};
                Master_tildeS2(s,d,thisLayer) = {eigenvalues};
                Master_Cprivate(s,d,thisLayer) = {privateCovarianceElements};
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute communication subspace
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [optDimReducedRankRegress,optLoss,optPredPerf_mean,optPredPerf_sem] = ...
                fct_compute_communicationSubspace(transpose(data_sender),transpose(data_receiver),...
                0:num_selNeur,numFolds,bool_print,bool_plot);
            
            % Write the data
            Master_ALLcomSub(s,d,:) = [optDimReducedRankRegress,optLoss,optPredPerf_mean,optPredPerf_sem];
            
            % Now that we know the optimal dimension of the communication subspace,
            % we perform Reduced-Rank Regression using this dimension alone
            [Brrr, B_, Gamma, Lambda] = ...
                ReducedRankRegress(transpose(data_receiver), transpose(data_sender), optDimReducedRankRegress);
            % B_ - Predictive dimensions: B_ = Bols*Gamma, where the columns of Gamma contain
            %       the eigenvectors of the optimal linear predictor Yhat = X*Bols (Bols is
            %       the ordinary least squared solution). The columns of B_, the
            %       predictive dimensions, are ordered according to target variance
            %        explained. For example, the top two predictive dimensions are B_(:,1:2).
            % Gamma  - The columns of Gamma contain the eigenvectors of the optimal linear
            %       predictor Yhat = X*Bols (Bols is the ordinary least squared solution).
            % Lambda - 	Associated eigenvalues of Yhat'*Yhat (the covariance
            %       matrix of Yhat).
            % Brrr - Reduced-rank regression solution to the problem: Y = XB
            
            % Write the data
            Master_FinalComSub(s,d,1) = {B_};
            Master_FinalComSub(s,d,2) = {Gamma};
            Master_FinalComSub(s,d,3) = {diag(Lambda)};
            Master_FinalComSub(s,d,4) = {Brrr};
            
        end
        
        % Save results
        thisCurrentParam = theseParams(p);
        filenameSave = strrep(sprintf([data_folder,'FAandComSub_50fromDisc_bin50ms_',thisParam,'%.03g_popRadius%.03g'],theseParams(p),pop_radius),'.','d');
        save(filenameSave,'thisParam','thisCurrentParam','paramName','paramUnit','colororder','colormapName',...
            'pop_radius','Master_pop_radius','pop_sizes','numDraws',...
            'Master_idxSelected','Master_meanFR',...
            'Master_sharedPR','Master_shared95',...
            'Master_Cshared','Master_W','Master_tildeU','Master_tildeS2','Master_Cprivate',...
            'Master_ALLcomSub','Master_FinalComSub')
    end
end
