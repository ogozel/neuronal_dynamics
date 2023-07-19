% Plot the prediction performance of the communication subspace as a
% function of the difference in PR in Sender and Receiver, as well as other
% related plots concerning interplay between E/I balance destabilization,
% size of disc radius from which neurons are sampled, prediction
% performance and dimensionality of the communication subspace.

% First run 'compute_FAcomSub.m' on cluster, save the results, and copy
% them locally. Then run this code on local machine.


clear all
% NB: comment following line to get figure(5000) with all the results for
% all different parameter sets on the same plot
close all 
clc


dataFolder = '../data_analysis/';

% Choose which parameter we are interested in
thisParam = 'tauIdecayinL1';
% 'sigmaRRfromIinL1'; 'sigmaRRfromIinL2';
% 'tauIdecayinL1'; 'tauIdecayinL2';

% Size of the disc from which neurons are sampled for the communication
% subspace prediction performance
% 1 is 5 neurons;   2 is 10 neurons;    3 is 20 neurons; 
% 4 is 30 neurons;  5 is 40 neurons;    6 is 100 neurons
discSize = 5; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(thisParam,'sigmaRRfromIinL1')
    
    theseParams = [0.1:0.025:0.3];
    paramName = '\sigma_{I}';
    paramUnit = '';
    colororder1 = winter(length(theseParams));
    colormapName = 'winter';
    markerShape = 'square';
    
elseif strcmp(thisParam,'sigmaRRfromIinL2')
    
    theseParams = [0.1:0.025:0.3];
    paramName = '\sigma_{I}';
    paramUnit = '';
    colororder1 = winter(length(theseParams));
    colormapName = 'winter';
    markerShape = 'diamond';
    
elseif strcmp(thisParam,'tauIdecayinL1')
    
    theseParams = [8:2:24];
    paramName = '\tau_{Id}';
    paramUnit = 'ms';
    colororder1 = copper(length(theseParams));
    colormapName = 'copper';
    markerShape = 'square';
    
elseif strcmp(thisParam,'tauIdecayinL2')
    
    theseParams = [8:2:24];
    paramName = '\tau_{Id}';
    paramUnit = 'ms';
    colororder1 = copper(length(theseParams));
    colormapName = 'copper';
    markerShape = 'diamond';
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover the within-area dimensionality and communication subspace
% measures computed with 'compute_FAcomSub.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CS_meanDimPerPopR = [];
CS_semDimPerPopR = [];
CS_meanPredPerfPerPopR = [];
CS_semPredPerfPerPopR = [];

meanPR_L1 = [];
meanPR_L2 = [];
semPR_L1 = [];
semPR_L2 = [];

for p=1:length(theseParams)
    
    filename = strrep(sprintf(['FAandComSub_50fromDisc_bin50ms_',...
        thisParam,'%.03g_popRadius100'],theseParams(p)),'.','d');
    
    % Load the daza
    dataS = load([dataFolder,filename,'.mat']);
    
    % Radius of the discs from which neurons were sampled
    MasterPopulationR = dataS.Master_pop_radius;
    
    % dataS.Master_ALLcomSub is of size: 
    %   (disc radiuses from which neurons are sampled)x(number of draws)x4
    % Slice 1: dimensionality of communication subspace
    % Slice 2: optimal loss for reduced-rank regression
    % Slice 3: mean prediction performance of communication subspace
    % Slice 4: SEM of prediction performance of communication subspace
    CS_AvgOverDraws = squeeze(mean(dataS.Master_ALLcomSub,2));
    CS_SEMOverDraws = squeeze(std(dataS.Master_ALLcomSub,1,2))/sqrt(size(dataS.Master_ALLcomSub,2));
    meanPR_AvgOverDraws = squeeze(mean(dataS.Master_sharedPR,2));
    semPR_AvgOverDraws = squeeze(std(dataS.Master_sharedPR,1,2))/sqrt(size(dataS.Master_sharedPR,2));
    
    % Combine data about communication subspace
    CS_meanDimPerPopR = cat(2, CS_meanDimPerPopR, CS_AvgOverDraws(:,1));
    CS_semDimPerPopR = cat(2, CS_semDimPerPopR, CS_SEMOverDraws(:,1));
    CS_meanPredPerfPerPopR = cat(2, CS_meanPredPerfPerPopR, CS_AvgOverDraws(:,3));
    CS_semPredPerfPerPopR = cat(2, CS_semPredPerfPerPopR, CS_SEMOverDraws(:,4));
    
    % Combine data about within-area dimensionality
    meanPR_L1 = cat(2, meanPR_L1, meanPR_AvgOverDraws(:,1));
    meanPR_L2 = cat(2, meanPR_L2, meanPR_AvgOverDraws(:,2));
    semPR_L1 = cat(2, semPR_L1, semPR_AvgOverDraws(:,1));
    semPR_L2 = cat(2, semPR_L2, semPR_AvgOverDraws(:,2));
    
    
    % Shared PR as a function of the radius of the disc population from
    % which neurons are sampled (in neurons)
    figure(101)
    hold on
    errorbar(MasterPopulationR/200,mean(dataS.Master_sharedPR(:,:,1),2),...
        std(dataS.Master_sharedPR(:,:,1),1,2)/sqrt(size(dataS.Master_sharedPR,2)),...
        'color',colororder1(p,:),'Linewidth',2,...
        'DisplayName',[paramName,' = ',num2str(theseParams(p)),paramUnit])
    xlabel('Disc radius','interpreter','latex')
    ylabel('Shared PR','interpreter','latex')
    title('L1','interpreter','latex')
    pbaspect([1 1 1])
    box off
    set(gca,'TickDir','out','Fontsize',18,'TickLabelInterpreter','latex')
    xticks(union(MasterPopulationR/200,0.25))
    yticks(0:5:25)
    ylim([0 25])
    
    figure(102)
    hold on
    errorbar(MasterPopulationR/200,mean(dataS.Master_sharedPR(:,:,2),2),...
        std(dataS.Master_sharedPR(:,:,2),1,2)/sqrt(size(dataS.Master_sharedPR,2)),...
        'color',colororder1(p,:),'Linewidth',2,...
        'DisplayName',[paramName,' = ',num2str(theseParams(p)),paramUnit])
    xlabel('Disc radius','interpreter','latex')
    ylabel('Shared PR','interpreter','latex')
    title('L2','interpreter','latex')
    pbaspect([1 1 1])
    box off
    set(gca,'TickDir','out','Fontsize',18,'TickLabelInterpreter','latex')
    xticks(union(MasterPopulationR/200,0.25))
    yticks(0:5:25)
    ylim([0 25])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction performance of the communication subspace as a function of
% disc radius from which neurons are sampled for the network with standard
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
k=1;
for r=1:length(MasterPopulationR)
    thisDiscR = MasterPopulationR(r);
    thisMean_predPerf = CS_meanPredPerfPerPopR(r,k);
    thisSEM_predPerf = CS_semPredPerfPerPopR(r,k);
    errorbar(thisDiscR,thisMean_predPerf,thisSEM_predPerf,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
end
xlabel('Population disc radius')
ylabel('Prediction performance [%] (\mu \pm SEM)')
title('Standard parameters')
box off
set(gca,'TickDir','out','FontSize',12)
colormap(colormapName)
c = colorbar('Ticks',linspace(0,1,length(theseParams)),'TickLabels',theseParams);
c.Label.String = paramName;
c.Label.FontSize = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction performance of the communication subspace as a function of
% parameter change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    thisSimParam = theseParams(k);
    thisMean_predPerf = CS_meanPredPerfPerPopR(discSize,k);
    thisSEM_predPerf = CS_semPredPerfPerPopR(discSize,k);
%     errorbar(thisSimParam,thisMean_predPerf,thisSEM_predPerf,...
%         'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
%         'color',colororder1(k,:),'Linewidth',2)
    plot(thisSimParam,thisMean_predPerf,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(ceil(length(theseParams)/2),:),'Linewidth',2,...
        'Marker',markerShape)
end
xlabel(thisParam)
ylabel('Prediction performance [%] (\mu \pm SEM)')
box off
set(gca,'TickDir','out','FontSize',12)
% colormap(colormapName)
% c = colorbar('Ticks',linspace(0,1,length(theseParams)),'TickLabels',theseParams);
% c.Label.String = paramName;
% c.Label.FontSize = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction performance as a function of PR in Sender (L1), or PR in
% Receiver (L2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    thisPRinL1 = meanPR_L1(discSize,k); % within-area PR in L1
    thisPRinL2 = meanPR_L2(discSize,k); % within-area PR in L2
    thisMean_predPerf = CS_meanPredPerfPerPopR(discSize,k);
    thisSEM_predPerf = CS_semPredPerfPerPopR(discSize,k);
    
    subplot(1,2,1)
    hold on
    errorbar(thisPRinL1,thisMean_predPerf,thisSEM_predPerf,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
    xlabel('Sender shared PR')
    ylabel('Prediction performance [%] (\mu \pm SEM)')
    box off
    
    subplot(1,2,2)
    hold on
    errorbar(thisPRinL2,thisMean_predPerf,thisSEM_predPerf,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
    xlabel('Receiver shared PR')
    ylabel('Prediction performance [%] (\mu \pm SEM)')
    box off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR in Sender (L1), or PR in Receiver (L2) as a function of the modified
% parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    thisSimParam = theseParams(k);
    thisPRinL1 = meanPR_L1(discSize,k); % within-area PR in L1
    thisPRinL2 = meanPR_L2(discSize,k); % within-area PR in L2
    thissemPRinL1 = semPR_L1(discSize,k);
    thissemPRinL2 = semPR_L2(discSize,k);
    thisMean_dimCS = CS_meanDimPerPopR(discSize,k);
    thisSEM_dimCS = CS_semDimPerPopR(discSize,k);
    
    errorbar(thisSimParam,thisPRinL1,thissemPRinL1,...
        'DisplayName','Sender shared PR',...
        'color','cyan','Linewidth',2)
    errorbar(thisSimParam,thisPRinL2,thissemPRinL2,...
        'DisplayName','Receiver shared PR',...
        'color','magenta','Linewidth',2)
    errorbar(thisSimParam,thisMean_dimCS,thisSEM_dimCS,...
        'DisplayName','Receiver shared PR',...
        'color','k','Linewidth',2)
end
xlabel(thisParam)
ylabel('Shared PR')
legend('Sender','Receiver','ComSub')
box off
set(gca,'TickDir','out','FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction performance as a function of PR in Sender (L1) minus PR in
% Receiver (L2)
% NB: if not 'close all', run for each parameter set to get all results on
% the same plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5000)
hold on
for k=1:length(theseParams)
    thisPRinL1 = meanPR_L1(discSize,k); % within-area PR in L1
    thisPRinL2 = meanPR_L2(discSize,k); % within-area PR in L2
    thisMean_predPerf = CS_meanPredPerfPerPopR(discSize,k);
    thisSEM_predPerf = CS_semPredPerfPerPopR(discSize,k);
    
%     errorbar(thisPRinL1-thisPRinL2,thisMean_predPerf,thisSEM_predPerf,...
%         'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
%         'color',colororder1(k,:),'Linewidth',2)
    plot(thisPRinL1-thisPRinL2,thisMean_predPerf,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(ceil(length(theseParams)/2),:),'Linewidth',2,...
        'Marker',markerShape)
end
xlabel('Sender shared PR - Receiver shared PR')
ylabel('Prediction performance [%] (\mu \pm SEM)')
box off
set(gca,'TickDir','out','FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensionality of communication subspace as a function of PR in Sender
% (L1) minus PR in Receiver (L2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    thisPRinL1 = meanPR_L1(discSize,k); % within-area PR in L1
    thisPRinL2 = meanPR_L2(discSize,k); % within-area PR in L2
    thisMean_dimCS = CS_meanDimPerPopR(discSize,k);
    thisSEM_dimCS = CS_semDimPerPopR(discSize,k);
    
    errorbar(thisPRinL1-thisPRinL2,thisMean_dimCS,thisSEM_dimCS,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
end
xlabel('Sender shared PR - Receiver shared PR')
ylabel('Dim comSub (\mu \pm SEM)')
box off
set(gca,'TickDir','out','FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensionality of communication subspace as a function of size of disc
% from which neurons are sampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    
    thisMean_dimCS = CS_meanDimPerPopR(:,k);
    thisSEM_dimCS = CS_semDimPerPopR(:,k);
    
    errorbar(MasterPopulationR,thisMean_dimCS,thisSEM_dimCS,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
end
xlabel('Population disc radius')
ylabel('Dim comSub (\mu \pm SEM)')
box off
set(gca,'TickDir','out','FontSize',12)
xticks(union(MasterPopulationR,50))
ylim([0 35])
title(thisParam)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensionality of communication subspace as a function of prediction
% performance of the communication subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for k=1:length(theseParams)
    
    thisMean_predPerf = CS_meanPredPerfPerPopR(:,k);
    thisMean_dimCS = CS_meanDimPerPopR(:,k);
    thisSEM_dimCS = CS_semDimPerPopR(:,k);
    
    errorbar(thisMean_predPerf,thisMean_dimCS,thisSEM_dimCS,...
        'DisplayName',[paramName,' = ',num2str(theseParams(k)),paramUnit],...
        'color',colororder1(k,:),'Linewidth',2)
end
xlabel('Mean predPerf of comSub')
ylabel('Dim comSub (\mu \pm SEM)')
box off
set(gca,'TickDir','out','FontSize',12)







