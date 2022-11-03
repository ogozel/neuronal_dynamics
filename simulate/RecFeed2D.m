% Code based on Chengcheng Huang's code:
% https://github.com/hcc11/SpatialNeuronNet
% https://github.com/hcc11/FI_SpatialNet
% I modified and drastically simplified the code for my project

function RecFeed2D(varargin)
% RecFeed2D(option, ParamChange)

% param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr,
%       gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
%       maxns, Irecord, Psyn
%   Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
%   Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii]; % out-degrees are fixed
%   taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
%   sigmaRR=[sigmaee, sigmaei; sigmaie, sigmaii];
%   sigmaRX=[sigmaeX; sigmaiX];

%   Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices,
%       sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
%       cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population.
%       For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and
%       Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop.
%   Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
%       The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.

%   conversion of neuron ID (exc) to (x,y) coordinate in [1, Ne1]x[1, Ne1]:
% exc. ID [1, Ne], x=ceil(I/Ne1); y=(mod((I-1),Ne1)+1); ID=(x-1)*Ne1+y
% inh. ID [Ne+1, Ne+Ni], x=ceil((I-Ne)/Ni1); y=(mod((I-Ne-1),Ni1)+1); ID=(x-1)*Ni1+y+Ne;

% sx: spike trains from Layer0
% s1: spike trains from Layer1 (only excitatory neurons)
% s2: spike trains from Layer2 (only excitatory neurons)
%     sg(1,:) contains spike times (for g=x,1,2)
%     sg(2,:) contains indices of neurons that spike (for g=x,1,2)
%
% data save in filename
% options is a struct w/ fields:
%   'save','CompCorr','fixW','savecurrent','loadS1'. Default values are 0.
% ParamChange is a cell of 2 columns,
%    the 1st column is variable names and
%    the 2nd column is the values.

% if options.save or options.savecurrent is 1, ParamChange needs to have field 'filename'.
% if options.CompCorr is 1, ParamChange needs to have field 'Nc',e.g. Nc=[500 500];
%      # of neurons to sample from Layer1 & Layer2.
% if options.fixW is 1, ParamChange needs to have field 'Wseed1' & 'Wseed2'.

nVarargs = length(varargin);
switch nVarargs
    case 1
        option = varargin{1};
    case 2
        option = varargin{1};
        ParamChange = varargin{2};
end

if ~isfield(option, 'save') option.save=0; end
if ~isfield(option, 'CompCorr') option.CompCorr=0; end
if ~isfield(option, 'loadS1') option.loadS1=0; end
if ~isfield(option, 'fixW') option.fixW=0; end
if ~isfield(option, 'useWfile') option.useWfile=0; end
if ~isfield(option, 'savecurrent') option.savecurrent=0; end
if ~isfield(option, 'saveRm') option.saveRm=0; end
if ~isfield(option, 'saveSx') option.saveSx=0; end
if ~isfield(option, 'saveS1') option.saveS1=1; end
if ~isfield(option, 'saveS2') option.saveS2=1; end
% Added by OG: save spike counts or not (E1, I1, E2, I2)
if ~isfield(option, 'saveSpkCounts') option.saveSpkCounts=0; end
if ~isfield(option, 'saveParam') option.saveParam=0; end
if ~isfield(option, 'saveW') option.saveW=0; end

if option.save==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end
if option.CompCorr==1
    if ~ismember('Nc',ParamChange(:,1))
        error('No Nc (1x2): # of neurons to sample to compute correlations')
    end
end
if option.loadS1==1
    if ~ismember('s1_fname',ParamChange(:,1))
        error('No s1_fname')
    end
end
if option.useWfile==1
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.saveW==1
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.savecurrent==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end

%% define parameters
dim ='2D';
% Number of neurons in network
Ne11=200; % Number in each direction
Ni11=100;
Ne21=200; % Number in each direction
Ni21=100;
Nx1=50; % input layer

param(1).Ne=Ne11*Ne11;
param(1).Ni=Ni11*Ni11;
param(1).Nx=Nx1*Nx1;
param(2).Ne=Ne21*Ne21;
param(2).Ni=Ni21*Ni21;
param(2).Nx=Ne11*Ne11;
param(1).N=param(1).Ne+param(1).Ni;
param(2).N=param(2).Ne+param(2).Ni;

% Stimulus param (by default, Poisson input)
p_stim.Nstim=1;
p_stim.stim_type='Uncorr';

% Rate of neurons in feedforward layer (kHz), size 1xNstim cell of element 
% 1xNsource array
p_stim.rX=.01;

% Number of sources for global correlation, size 1xNstim
p_stim.Nsource=1;  

% Total simulation time (in msec)
T=20000;
% Bin size in [ms], timestep for the forward Euler method
dt=.05;  
Tburn=1000; % Burn-in period

% Static currents to Layer 2 (none here)
inE=0;
inI=0;

% Connection widths
param(1).sigmaRX=.05*ones(2,1);
param(1).sigmaRR=.1*ones(2,2);
param(2).sigmaRX=.05*ones(2,1);
param(2).sigmaRR=.1*ones(2,2);

% Number of neurons to record synaptic inputs and voltages from
nrecordE0=zeros(1,2);
nrecordI0=zeros(1,2);

% Synaptic time constants
% 1st row: feedforward synapses from L1 to L2
% 2nd row: recurrent synapses from L2 excitatory neurons
% 3rd row: recurrent synapses from L2 inhibitory neurons
param(2).taudsyn=[5; 5; 8];
param(2).taursyn=[1; 1; 1];
param(2).Psyn=[1; 1; 1];

% Same parameters for L1
param(1).taudsyn=[5; 5; 8];
param(1).taursyn=[1; 1; 1];
param(1).Psyn=[1; 1; 1];

% Mean connection probabilities
% Prr = [Pee Pei; Pie, Pii], Prx = [Pex; Pix]
param(1).Prr=[.01, .04; .03, .04];
param(2).Prr=[.01, .04; .03, .04];
param(1).Prx=[ .1; .05]; % from L0 (input) to L1
param(2).Prx=[ .1; .05]; % from L1 to L2

% Connection strengths (scaled by sqrt(N) later)
% Recurrent conn strengths in [mV], Jr=[Jee Jei; Jie, Jii]
param(1).Jr=[80 -240; 40, -300];
param(2).Jr=[80 -240; 40, -300];
% Feedforward connection strengths to L1 and L2 in [mV], Jx=[JeF; JiF]
param(1).Jx=[240; 400]; % from L0 (input layer) to L1
param(2).Jx=[15; 25]; % from L1 to L2
% OG: values modified to be equivalent to L1 (16 times more feedforward 
% input = 40'000/2'500 exc. neur)

param(1).Iapp = cat(1,0*ones(param(1).Ne,1),0*ones(param(1).Ni,1));
param(2).Iapp = cat(1,inE*ones(param(2).Ne,1),inI*ones(param(2).Ni,1));

% Change parameters
if nVarargs==2
    for i=1:size(ParamChange,1)
        eval([ParamChange{i,1} '= ParamChange{i,2};']);
    end
end

if option.loadS1
    p_stim.s1_fname=s1_fname;
end

clear ParamChange varargin;

%% Initialization

% Max number of spikes is 50'000'000 in L1 and L2
maxrate=[.05 .05]; 
param(1).maxns=param(1).N*T*maxrate(1);
param(2).maxns=param(2).N*T*maxrate(2);
fprintf('\nmaximum average rates to record: layer 1 %d, layer 2 %d\n (Hz)',round(maxrate(1)*1e3), round(maxrate(2)*1e3))

for par=1:2
    param(par).dt=dt;
    param(par).T=T;
    % EIF neuron paramters
    param(par).gl=[1/15 1/10];  % E, I
    param(par).Cm=[1 1];
    param(par).vlb=[-100 -100];
    param(par).vth=[-10 -10];
    param(par).DeltaT=[2 .5];
    param(par).vT=[-50 -50]; %mV
    param(par).vre=[-65 -65];
    param(par).tref=[1.5 .5];
    V0min=param(par).vre(1);
    V0max=param(par).vT(1);
    % OG: modified so it can be different for each neuron
    param(par).vl=param(par).Iapp.*cat(1, 15*ones(param(par).Ne,1),...
        10*ones(param(par).Ni,1))-60*ones(param(par).N,1);
    param(par).V0=(V0max-V0min).*rand(param(par).N,1)+V0min;
    % Same recurrent out-degrees for L1 and L2
    % Kee=400, Kei=1600, Kie=300, Kii=400
    param(par).Kr=ceil(param(par).Prr.*[param(par).Ne,...
        param(par).Ne; param(par).Ni,param(par).Ni]);
    % Different feedforward out-degrees for L1 and L2
    % From L0 (input) to L1: KeF=4000, KiF=500
    % From L1 to L2: KeF=2000, KiF=500
    param(par).Kx=ceil(param(par).Prx.*[param(par).Ne; param(par).Ni]);
    % Neuron indice to record synaptic currents and Vm
    param(par).Irecord=[randi(param(par).Ne,1,nrecordE0(par)),...
        (randi(param(par).Ni,1,nrecordI0(par))+param(par).Ne)];
    param(par).Jr=param(par).Jr/sqrt(param(par).N);
    param(par).Jx=param(par).Jx/sqrt(param(par).N);
    
    % Effective connection weights
    q=param(par).Ne/param(par).N;
    wrx=(param(par).Jx(:,1)).*param(par).Prx*param(par).Nx/param(par).N;
    wrr=(param(par).Jr).*param(par).Prr.*[q, 1-q; q, 1-q];
    % For balanced state to exist this vector should be decreasing
    fprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n\n',...
        wrx(1)/wrx(2),abs(wrr(1,2)/wrr(2,2)),abs(wrr(1,1)/wrr(2,1)));
    % and these values should be >1
    fprintf('\nAlso, this number should be greater than 1: %.2f\n\n',...
        abs(wrr(2,2)/wrr(1,1)));
end
scurr = rng;

%% generate input spike trains
if option.loadS1==0
    sx=genXspk(p_stim,param(1).Nx,T);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%         Simulation              %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((param(2).N)<=200000)
    disp('simulation starts')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set-up Layer 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('set-up Layer1')
    % Random initial membrane potentials
    if option.loadS1
        disp('load s1')
        load(s1_fname,'s1');
        disp(['In the case that option.loadS1=1, the following line ',...
            'should avoid ID=0 in s1(2,:), and saves only the spikes ',...
            'of excitatory neurons!'])
        % Only the spikes of excitatory neurons are saved!
        s1=s1(:,s1(2,:)<=param(1).Ne+.1&s1(2,:)>0&s1(1,:)<T);
        if option.save
            save(filename,'T')
        end
    else
        disp('We are in the case loadS1=0.')
        if option.fixW
            disp('We are in the case fixW=1.')
            if option.useWfile==1
                load(W_fname,'Wrr1','Wrf1')
                fprintf('load weight from %s\n',W_fname)
            else
                param(1).Wseed=Wseed1;
                rng(Wseed1,'twister')
                fprintf('seed%d for Wrr1, Wrf1\n',Wseed1)
                disp('generating :q!, Wrf1')
                tic
                [Wrr1,Wrf1]=gen_weights(param(1).Ne,param(1).Ni,...
                    param(1).Nx,param(1).sigmaRX,param(1).sigmaRR,...
                    param(1).Prr, param(1).Prx,dim);
                elapsetime=toc;
                fprintf('elapsetime=%.2f sec\n',elapsetime)
            end
        else
            disp('We are in the case fixW=0.')
            disp('generating Wrr1, Wrf1')
            tic
            Wseed1 = rng;
            [Wrr1,Wrf1]=gen_weights(param(1).Ne,param(1).Ni,param(1).Nx,...
                param(1).sigmaRX,param(1).sigmaRR,param(1).Prr,...
                param(1).Prx,dim);
            elapsetime=toc;
            fprintf('elapsetime=%.2f sec\n',elapsetime)
        end
        if option.saveW
            disp('We have saveW=1, so we save W_fname for L1.')
            save(W_fname,'Wrr1','Wrf1','Wseed1','param')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulate Layer 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('simulating Layer1')
        tic
        [s1,Isyn1,Vm1]=EIF1DRFfastslowSyn(sx,Wrf1,Wrr1,param(1));
        clear Wrr1 Wrf1;
        
        End=find(s1(2,:)==0,1)-1;
        disp(['In the case that option.loadS1=0, the following line ',...
            'should avoid ID=0 in s1(2,:) BEFORE s1 is saved!'])
        s1 = s1(:,s1(2,:)~=0);
        elapsetime=toc;
        nuSim(3)=1000*nnz(s1(1,1:End)>Tburn & s1(2,1:End)<=param(1).Ne)/(param(1).Ne*(T-Tburn));  % Hz
        nuSim(4)=1000*nnz(s1(1,1:End)>Tburn & s1(2,1:End)>param(1).Ne)/(param(1).Ni*(T-Tburn));
        nuSim(5)=1000*nnz(sx(1,:)>Tburn & sx(1,:)<T)/(param(1).Nx*(T-Tburn));
        fprintf('\naverage rates \n E1: %.2f, I1: %.2f, X: %.2f \n elapsetime=%.2f sec\n',...
            nuSim(3),nuSim(4),nuSim(5),elapsetime)
        
        if option.save
            if option.saveSpkCounts % Save spike counts
                Tw=50; % timewindow in [ms] for the spike count computations
                Nt=floor(T/Tw); % number of non-overlapping timewindows
                idx=1:Nt;
                X=spktime2count(sx,1:param(1).Nx,Tw,Nt,1);
                if strcmp(p_stim.stim_type,'OriMap_gabor_Tseg')
                    Nt_on=p_stim.T_on/Tw;
                    Nt_off=p_stim.T_off/Tw;
                    Nseg=Nt_on+Nt_off;
                    idx(mod(idx-1,Nseg)+1<=Nt_off)=0;
                    Xstim=permute(sum(reshape(X(:,idx>.5),param(1).Nx,...
                        Nt_on,[]),2),[1 3 2]);
                else % if the stimulus is always on, keep Xstim as it is
                    Xstim=X;
                end
                save(filename,'T','Xstim','Tw')
                if option.saveSx
                    save(filename,'sx','-append')
                end
                clear Xstim X;
                E1=int8(spktime2count(s1(:,1:End),1:param(1).Ne,Tw,Nt,1));
                I1=int8(spktime2count(s1(:,1:End),...
                    (1+param(1).Ne):param(1).N,Tw,Nt,1));
                save(filename,'E1','I1', '-append')
                clear E1 I1
            else
                save(filename,'T')
                if option.saveSx
                    save(filename,'sx','-append')
                end
                
            end
            if option.saveRm
                re1=hist(s1(1,s1(2,:)<=param(1).Ne&s1(2,:)>0),1:T)/param(1).Ne*1e3;
                ri1=hist(s1(1,s1(2,:)>param(1).Ne),1:T)/param(1).Ni*1e3;
                save(filename,'re1','ri1','-append')
                clear re1 ri1;
            end
            if option.saveS1
                if exist('s1_fname','var')
                    save(s1_fname,'s1')
                    if strcmp(p_stim.stim_type,'OriMap_gabor_Tseg')
                        th_id=p_stim.th_id;
                        save(s1_fname,'th_id','-append')
                    end
                else
                    save(filename,'s1','-append')
                end
            end
        end
        disp(['The following line makes sure that s1 is in the proper ',...
            'form to set as input of L2 (only excitatory neurons in s1)'])
        s1=s1(:,s1(2,:)<=param(1).Ne+.1);
    end % end of case loadS1==0
    
    
    if option.saveParam
        save(filename,'param','p_stim','-append')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set-up Layer 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('set-up Layer2')
    
    if option.fixW
        if option.useWfile==1
            load(W_fname,'Wrr2','Wrf2')
            fprintf('load weight from %s\n',W_fname)
        else
            param(2).Wseed=Wseed2;
            rng(Wseed2,'twister')
            fprintf('seed%d for Wrr2, Wrf2\n',Wseed2)
            disp('generating Wrr2, Wrf2')
            tic
            [Wrr2,Wrf2]=gen_weights(param(2).Ne,param(2).Ni,param(2).Nx,...
                param(2).sigmaRX,param(2).sigmaRR,param(2).Prr,...
                param(2).Prx,dim);
            elapsetime=toc;
            fprintf('elapsetime=%.2f sec\n',elapsetime)
        end
    else
        disp('generating Wrr2, Wrf2')
        tic
        Wseed2 = rng;
        [Wrr2,Wrf2]=gen_weights(param(2).Ne,param(2).Ni,param(2).Nx,...
            param(2).sigmaRX,param(2).sigmaRR,param(2).Prr,...
            param(2).Prx,dim);
        elapsetime=toc;
        fprintf('elapsetime=%.2f sec\n',elapsetime)
    end
    if option.saveW
        disp('We have saveW=1, so we append L2 values to W_fname.')
        save(W_fname,'Wrr2','Wrf2','Wseed2','param','-append')
    end
    
    rng(scurr);
    Jx=param(2).Jx;
    
    % OG: modified so it can be different for each neuron
    param(2).vl=param(2).Iapp.*cat(1, 15*ones(param(2).Ne,1),...
        10*ones(param(2).Ni,1))-60*ones(param(2).N,1);
    param(2).V0=(V0max-V0min).*rand(param(2).N,1)+V0min;
    param(2).Jx=Jx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Layer 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('simulating Layer2')
    tic
    [s2,Isyn2,Vm2]=EIF1DRFfastslowSyn(s1, Wrf2,Wrr2,param(2));
    elapsetime=toc;
    
    End=find(s2(2,:)==0,1)-1;
    
    if option.saveRm
        re2=hist(s2(1,s2(2,:)<=param(2).Ne&s2(2,:)>0),1:T)/param(2).Ne*1e3;
        ri2=hist(s2(1,s2(2,:)>param(2).Ne),1:T)/param(2).Ni*1e3;
        save(filename,'re2','ri2','-append')
    end
    nuSim(1)=1000*nnz(s2(1,1:End)>Tburn & s2(2,1:End)<=param(2).Ne)/(param(2).Ne*(T-Tburn));  % Hz
    nuSim(2)=1000*nnz(s2(1,1:End)>Tburn & s2(2,1:End)>param(2).Ne)/(param(2).Ni*(T-Tburn));
    
    fprintf('\naverage rates \n E2: %.2f, I2: %.2f,\n elapsetime=%.2f sec\n',...
        nuSim(1),nuSim(2),elapsetime)
    
    if option.save
        if option.saveSpkCounts % save spike counts if needed
            Tw=50;
            Nt=floor(T/Tw);
            E2=int8(spktime2count(s2(:,1:End),1:param(2).Ne,Tw,Nt,1));
            I2=int8(spktime2count(s2(:,1:End),...
                (1+param(2).Ne):param(2).N,Tw,Nt,1));
            save(filename,'E2','I2','-append')
        end
        if option.saveS2
            s2=s2(:,s2(2,:)~=0);
            save(filename,'s2','nuSim','-append')
        end
    end
    
    param(2).Jx=Jx;
    
    clear Wrr2 Wrf2;
else
    error('N too large') % Your computer probably can't handle this
end
disp('simulation ends')

%% compute nearby covariance

if option.CompCorr==1
    rng('shuffle');
    [C,COV_d,Cbar,COVbar,daxis,rate1,rate2,var1,var2]=...
        corr_d(s1,s2,param(2).Nx,param(2).Ne,dim,Nc);
    fprintf('\nCee=%.4f, Cex=%.4f,Cxx=%.4f\n\n', C(1,1),C(1,2),C(1,3))
    Cbar
    %     FF=mean(var2./rate2)/5
    FF = [mean(var1./rate1)/5, mean(var2./rate2)/5]
    if option.save
        save(filename,'C','COV_d','Cbar','COVbar','daxis','rate1',...
            'rate2','var1','var2','FF','-append') % also save Fano factor
    end
end

if option.savecurrent
    Nrecord1=nrecordE0(1)+nrecordI0(1);
    Nskip=10; % record every 10 time steps
    IsynX1=param(1).Psyn(1)*Isyn1(1:Nrecord1,Nskip:Nskip:end);
    IsynE1=param(1).Psyn(2)*Isyn1((Nrecord1+1):2*Nrecord1,Nskip:Nskip:end);
    IsynI1=param(1).Psyn(3)*Isyn1((2*Nrecord1+1):3*Nrecord1,Nskip:Nskip:end);
    Vm1=Vm1(:,Nskip:Nskip:end);
    save(filename,'IsynX1','IsynE1','IsynI1','Vm1','-append')
    
    param(2).Psyn=[.2 .8; 1, 0; 1, 0];
    syntype_id=find(param(2).Psyn(:)>0);
    % type 1: X, 2:E, 3:I, for updating postsyn input type *
    syntype=mod(syntype_id-1,3)+1; 
    
    Nrecord2=nrecordE0(2)+nrecordI0(2);
    Nt=floor(size(Isyn2,2)/Nskip);
    IsynX2=zeros(Nrecord2,Nt);
    IsynE2=zeros(Nrecord2,Nt);
    IsynI2=zeros(Nrecord2,Nt);
    for ntype=1:length(syntype)
        switch  syntype(ntype)
            case 1
                IsynX2=IsynX2+...
                    param(2).Psyn(syntype_id(ntype))*...
                    Isyn2(((ntype-1)*Nrecord2+1):ntype*Nrecord2,Nskip:Nskip:end);
            case 2
                IsynE2=IsynE2+...
                    param(2).Psyn(syntype_id(ntype))*...
                    Isyn2(((ntype-1)*Nrecord2+1):ntype*Nrecord2,Nskip:Nskip:end);
            case 3
                IsynI2=IsynI2+...
                    param(2).Psyn(syntype_id(ntype))*...
                    Isyn2(((ntype-1)*Nrecord2+1):ntype*Nrecord2,Nskip:Nskip:end);
        end
    end
    Vm2=Vm2(:,Nskip:Nskip:end);
    save(filename,'IsynX2','IsynE2','IsynI2','Vm2','-append')
end

