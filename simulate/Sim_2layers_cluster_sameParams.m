% Code based on Chenhgcheng Huang's code:
% https://github.com/hcc11/SpatialNeuronNet
% https://github.com/hcc11/FI_SpatialNet
% I did some merging of the codes, as well as modifications in order to
% have the type of stimulations and outputs needed for my project

% This code runs simulations of a 3-layer spiking recurrent neural network
% L0: input layer
% L1: Sender network
% L2: Receiver network
% First of all, we need to compile the c code first:
% mex EIF1DRFfastslowSyn.c
% mex spktime2count.c



clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool_cluster = 1;

% Stimulation type
% 'Uncorr': Poisson input
% 'OriMap_gabor_Tseg' (with parameters specified below): for Gabor input
this_stimulationType = 'Uncorr';

% Network type
this_networkType = 'RecFeed2D';

% Measures we compute
bool_fixW = 1;
% if fixW=0 we get random weights; 
% if fixW=1 ParamChange needs Wseed1 and Wseed2
bool_compCorr = 1; % if 1, compute correlations (by default it is 0)
bool_saveSpkCnts = 1; % if 1, save spike counts from all layers
% (default is 0, do not save)

% Simulations
T = 20000;  % total time of simulation (ms)

% There are a total of x=(Ntrial x Np) trials
% In bash script (runslurm), it has to match: #SBATCH -a 1-x%y
% y is the number of jobs that are run in parallel
Ntrial = 30; % number of trials per parameter set
Np = 9; %18; % number of parameter sets


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is for setting parameters differently if running as a job array
% on cluster or as a simple job on local machine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bool_cluster % if on cluster
    
    rng('shuffle');
    AI = getenv('SLURM_ARRAY_TASK_ID');
    job_dex = str2num(AI) % number of the simulation
    % job_dex=1 to Ntrial correpond to sims with the first parameter set,
    % job_dex=(Ntrial+1) to (2*Ntrial) correspond to sims with the second
    % parameter set
    seed_offset = randi(floor(intmax/10));
    rng(job_dex + seed_offset);
    
    data_folder='/user_data/ogozel/';
    
else % if on local machine
    
    job_dex = 1;
    seed_offset = 1;
    
    data_folder='../data_sim/';
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial number (only goes from 1 to Ntrial)
nt = mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1);

% Parameter set for the current trial number
pid = mod(ceil(job_dex/Ntrial)-1,Np)+1;

Wseed1_range=[8541; 45134; 48395; 3547; 14845;  71109; 99911; 98570;...
    68790; 16203 ];

Wseed2_range=[800281; 141887; 421762; 915736; 792208; 959493; ...
    157614;  970593; 957167; 485376];
nws=3;


% Make sure that the orientation map exists, and if not create and save it
if exist([data_folder 'theta_map.mat'],'file')~=0
    load([data_folder 'theta_map.mat'],'theta_map')
else
    Nx=50;
    Kc=5*2*pi; % spatial frequency
    theta_map=ori_map(Nx,Kc);  % generate orientation map (size Nx by Nx)
    save([data_folder 'theta_map.mat'],'theta_map')
end

nType = pid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set stimulus type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(this_stimulationType,'OriMap_gabor_Tseg')
    
    disp('Gabor stimulation')
    
    testp.theta0 = 0.75; % only 1 orientation!!
    %     dtheta = 0.01;
    %     testp.theta0 = [.5-dtheta/2, .5+dtheta/2]; % orientations to test
    % testp.theta0=0.02:.02:1;
    Nth=length(testp.theta0); % Number of Gabor orientations
    
    % L4 neurons (input neurons, L0)
    rX=.01; % mean input rate in [kHz]
    dx=0.04;
    x=-.48:dx:.48;
    sigma=0.2;
    lambda=.6;
    % Noise
    sigma_n=3.5; % input noise intensity
    tau_n=40;  % tau of noise
    
elseif strcmp(this_stimulationType,'Uncorr')
    disp('Poisson stimulation')
else
    disp('Neither a Gabor nor Poisson stimulation')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set network type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Modify sigma_I_rec %%%
if pid<=9
    
    disp(['Modify sigmaRRfromI in L1 (Sender) to ',...
        num2str(sigmaRRfromI_range(nType))])
    
    nType=pid;
    sigmaRRfromI_range = [0.1:0.025:0.3];
    
    % Modification in L1
    param(1).sigmaRR(1,2) = sigmaRRfromI_range(nType);
    param(1).sigmaRR(2,2) = sigmaRRfromI_range(nType);
    
    Type=strrep(sprintf('L1_L2_sameParams_L1_sigmaRRFromI%.03g',...
        sigmaRRfromI_range(nType)),'.','d')
    
    opt.useWfile=1;
    this_W = ['weights_',Type];
    W_fname = sprintf('%s%c',data_folder,this_W)
    
else
    
    disp(['Modify sigmaRRfromI in L2 (Receiver) to ',...
        num2str(sigmaRRfromI_range(nType))])
    
    nType=pid-9;
    sigmaRRfromI_range = [0.1:0.025:0.3];
    
    % Modification in L2
    param(2).sigmaRR(1,2) = sigmaRRfromI_range(nType);
    param(2).sigmaRR(2,2) = sigmaRRfromI_range(nType);
    
    Type=strrep(sprintf('L1_L2_sameParams_L2_sigmaRRFromI%.03g',...
        sigmaRRfromI_range(nType)),'.','d')
    
    opt.useWfile=1;
    this_W = ['weights_',Type];
    W_fname = sprintf('%s%c',data_folder,this_W)
    
end

%%% Modify tau_I_decay %%%
% if pid<=9
%     
%     disp(['Modify tauIdecay in L1 (Sender) to ',...
%         num2str(tauIdecay_range(nType))])
%     
%     nType=pid;
%     tauIdecay_range = [8:2:24];
%     
%     % Modification in L1
%     param(1).taudsyn(3) = tauIdecay_range(nType);
%     
%     Type=strrep(sprintf('L1_L2_sameParams_L1_tauIdecay%.03g',...
%         tauIdecay_range(nType)),'.','d')
%     
%     % When modifying tauIdecay, we use the standard connectivity matrices 
%     opt.useWfile=1;
%     this_W = 'weights_L1_L2_sameParams';
%     W_fname = sprintf('%s%c',data_folder,this_W)
%       
% else
%     
%     nType=pid-9;
%     tauIdecay_range = [8:2:24];
%     
%     % Modification in L2
%     param(2).taudsyn(3) = tauIdecay_range(nType);
%     
%     Type=strrep(sprintf('L1_L2_sameParams_L2_tauIdecay%.03g',...
%         tauIdecay_range(nType)),'.','d')
%     
%     % When modifying tauIdecay, we use the standard connectivity matrices 
%     opt.useWfile=1;
%     this_W = 'weights_L1_L2_sameParams';
%     W_fname = sprintf('%s%c',data_folder,this_W)
%     
%     disp(['Modify tauIdecay in L2 (Receiver) to ',num2str(tauIdecay_range(nType))])
%     
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On the cluster, the different trials are run in parallel, so no actual
% need of the for loop, for loop only needed to test on local machine
for trial = nt
    rng(trial + seed_offset);
    
    % Decide which data to save
    opt.save=1; % save data
    opt.saveSx=0; % save spike times from Layer 0 (input layer) (default is 0, save to filename)
    opt.saveS1=1; % save spike times from Layer 1 (default is 1, save to filename, or s1_fname if specified)
    opt.saveS2=1; % save spike times from Layer 2 (default is 1, save to filename)
    opt.saveSpkCounts = bool_saveSpkCnts;
    
    opt.CompCorr = bool_compCorr;
    Nc = [500 500]; % number of excitatory neurons sampled from Layers 1 and 2
    
    opt.loadS1=0;
    
    opt.fixW = bool_fixW;
    Wseed1=Wseed1_range(nws);
    Wseed2=Wseed2_range(nws);
    opt.saveW=0; % NB: good to save the weights if they are not fixed
    opt.savecurrent=0; % To save the synaptic input to some neurons (hard-coded in RecFeed2D.m)
    opt.saveRm=0;
    
    if trial==1
        opt.saveParam=1;
    else
        opt.saveParam=0;
    end
    
    if exist([data_folder strcat(this_W,'.mat')],'file')~=0
        disp('The weights file this_W already exists (but might just be in the workspace/not saved), so we set useWfile=1.')
        opt.useWfile=1; % use the one that was already created and saved
    else
        disp('The weights file this_W does not exist yet, so we set useWfile=0.')
        opt.useWfile=0; % create and save weight matrices to W_fname
        disp('saveW=1 if the W_fname does not exist yet.')
        opt.saveW = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stimulation input (firing rate of L0 neurons, the input layer)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p_stim.stim_type = this_stimulationType;
    
    if strcmp(this_stimulationType,'OriMap_gabor_Tseg')
        [X,Y]=meshgrid(x,x);
        X=X(:);Y=Y(:);
        % Gabor image (size 625x1)
        Imag=@(theta) (exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta)))...
            .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))));
        NI=numel(Imag(0));
        % Gabor filters for each L4 neuron (size(Nx^2 x NI)
        Filter=@(theta) exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta))...
            .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))))...
            /(sum(Imag(0).^2));
        % Gabor filters for each
        F=Filter(theta_map(:)'*pi)';
        fr=F*Imag(mean(testp.theta0)*pi);
        F=F/mean(fr)*rX;
        
        % Firing rate for each theta (size Nx^2 x Nth)
        fr=F*Imag(testp.theta0*pi);
        
        p_stim.F=F; % firing rate (Nx by Nx)
        p_stim.theta_map=theta_map;
        p_stim.NI=NI;
        p_stim.rX=rX;
        p_stim.fr=fr;
        p_stim.sigma_n=sigma_n;
        p_stim.sigma=sigma;
        p_stim.lambda=lambda;
        p_stim.T=T;
        p_stim.tau_n=tau_n;
        
        % Continuous input in Layer 0
        p_stim.T_on = T;
        p_stim.T_off = 0;
        p_stim.rX_off=.005; % firing rate during OFF intervals (kHz)
        % Number of stim presentations per trial in case the stimulus is 
        % ON/OFF). If it is always ON: 1.
        Nseg=ceil(T/(p_stim.T_on+p_stim.T_off));
        
        % Randomly select orientation id for each stim presentation
        p_stim.th_id=int8(randsample(Nth,Nseg,1));
        p_stim.theta0=testp.theta0;
        
    end
    
    % Filename for the saved simulation results
    filename=strrep(sprintf(strcat('%s',this_networkType,'_',this_stimulationType,'_saveSpkCnts',num2str(bool_saveSpkCnts),'_fixW',num2str(bool_fixW),'_%s_ID%.0f'),...
        data_folder,Type,trial),'.','d');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Keep track of the parameters that have been changed from their
    % default values in the structure ParamChange
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Modify sigma_I_rec %%%
    if pid<=9
        ParamChange={'filename',filename;'T',T;'Nc',Nc;...
            'param(1).sigmaRR(1,2)',param(1).sigmaRR(1,2);...
            'param(1).sigmaRR(2,2)',param(1).sigmaRR(2,2)};
    else
        ParamChange={'filename',filename;'T',T;'Nc',Nc;...
            'param(2).sigmaRR(1,2)',param(2).sigmaRR(1,2);...
            'param(2).sigmaRR(2,2)',param(2).sigmaRR(2,2)};
    end
    
    %%% Modify tau_I_decay %%%
%     if pid<=9
%         ParamChange={'filename',filename;'T',T;'Nc',Nc;...
%             'param(1).taudsyn(3)',param(1).taudsyn(3)};
%     else
%         ParamChange={'filename',filename;'T',T;'Nc',Nc;...
%             'param(2).taudsyn(3)',param(2).taudsyn(3)};
%     end
%     
    
    if strcmp(this_stimulationType,'OriMap_gabor_Tseg')
        ParamChange=[ParamChange;{'p_stim',p_stim}];
    end
    
    if opt.loadS1
        ParamChange=[ParamChange;{'s1_fname',s1_fname}];
    end
    
    if opt.fixW && opt.useWfile==0
        disp('fixW=1 and useWfile=0, so we write Wseed1 and Wseed2 in ParamChange.')
        ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
    end
    
    if opt.useWfile==1
        disp('OG: useWfile=1, so we check if it exists already.')
        if exist([W_fname '.mat'], 'file')==0
            disp('If W_fname does not exist already, saveW is set to 1.')
            opt.saveW=1;
        end
    end
    
    if opt.saveW==1||opt.useWfile==1
        disp('OG: if saveW=1 or useWfile=1, write W_fname in ParamChange.')
        ParamChange=[ParamChange;{'W_fname',W_fname}];
    end
    
    clear F Filter theta_map X Y Imag;
    
    disp('Write down filename:')
    filename
    RecFeed2D(opt, ParamChange) % main simulation
    
    if strcmp(this_stimulationType,'OriMap_gabor_Tseg')
        th_id=p_stim.th_id;
        save(filename,'testp','th_id','ParamChange','-append')
    else
        save(filename,'ParamChange','-append')
    end
    
end

