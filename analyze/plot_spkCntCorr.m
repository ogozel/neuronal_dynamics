% Plot the distribution of pairwise spike count correlations within L1, 
% within L2, and between L1 and L2, as well as the average pairwise 
% correlation as a function of pairwise distance.

% NB: If this script is run on a local machine, then one trial of the
% simulations (here the one with ID1) for each parameter set of interest 
% has to be copied on the local machine. Modify the path to the data
% accordingly.

clear all
close all
clc

% Add path to corr_d function and the data of interest
addpath('../simulate/')
addpath('../data_cluster_sigmaI/')
addpath('../data_cluster_tauIdecay/')


layerOfChange = 1; % 1 or 2
paramVal = 0.3; % 0.3; 0.1
filename = strrep(sprintf('RecFeed2D_Uncorr_saveSpkCnts1_fixW1_L1_L2_sameParams_L%.03g_sigmaRRFromI%.03g_spkCntCorr_ID1',...
        layerOfChange,paramVal),'.','d');
load(['../data_cluster_sigmaI/',filename,'.mat'])


N1 = length(rate1); % number of selected neurons in L1
N2 = length(rate2); % number of selected neurons in L2

corrL2 = Corr(1:N2,1:N2); % Layer2 is first
corrL1 = Corr(N2+1:end,N2+1:end);
corrL1L2 = Corr(1:N2,N2+1:end); % between-area correlation

corrL1 = triu(corrL1,1);
corrL1 = reshape(corrL1,[N1*N1,1]);
idx = find(corrL1~=0);
corrL1 = corrL1(idx);

corrL2 = triu(corrL2,1);
corrL2 = reshape(corrL2,[N2*N2,1]);
idx = find(corrL2~=0);
corrL2 = corrL2(idx);

corrL1L2 = reshape(corrL1L2,[N1*N2,1]);


%% Plotting

edges = -1:0.01:1;

figure()
[bincountsL1,edges] = histcounts(corrL1,edges);
histogram('BinCounts', bincountsL1/sum(bincountsL1), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['Layer 1 - sigmaI=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL1))])
box off
set(gca,'TickDir','out')

figure()
[bincountsL2,edges] = histcounts(corrL2,edges);
histogram('BinCounts', bincountsL2/sum(bincountsL2), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['Layer 2 - sigmaI=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL2))])
box off
set(gca,'TickDir','out')

figure()
[bincountsL1L2,edges] = histcounts(corrL1L2,edges);
histogram('BinCounts', bincountsL1L2/sum(bincountsL1L2), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['L1-L2 - sigmaI=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL1L2))])
box off
set(gca,'TickDir','out')

% Compute and plot the pairwise noise correlation as a function of pairwise
% distance
[C_d,COV_d,Cbar,COVbar,daxis,rate1,rate2,var1,var2]=corr_d(s1,s2,40000,40000,'2D',[500, 500]);
% C_d(:,1) : L2-L2
% C_d(:,2) : L1-L2
% C_d(:,3) : L1-L1

figure()
plot(daxis, C_d(:,3), 'k', 'Linewidth', 2)
hold on
yline(0)
xlabel('Pairwise distance')
ylabel('Spike count correlation')
box off
set(gca,'TickDir','out')
ylim([-0.2 1])


%% For temporal destabilization

layerOfChange = 2; % 1 or 2
paramVal = 24;
filename = strrep(sprintf('RecFeed2D_Uncorr_saveSpkCnts1_fixW1_L1_L2_sameParams_L%.03g_tauIdecay%.03g_ID1',...
        layerOfChange,paramVal),'.','d');
load(['../data_cluster_tauIdecay/',filename,'.mat'])

[C_d,COV_d,Cbar,COVbar,daxis,rate1,rate2,var1,var2]=corr_d(s1,s2,40000,40000,'2D',[500, 500]);

figure()
plot(daxis, C_d(:,3), 'k', 'Linewidth', 2)
hold on
yline(0)
xlabel('Pairwise distance')
ylabel('Spike count correlation')
box off
set(gca,'TickDir','out')
ylim([-0.2 1])

% Compute the pairwise correlations
[Corr, Cov, rate2, rate1, var2, var1] = fct_spkCntCorr(s1,s2,40000,40000,'2D',[500, 500]);

N1 = length(rate1); % number of selected neurons in L1
N2 = length(rate2); % number of selected neurons in L2

corrL2 = Corr(1:N2,1:N2); % Layer2 is first
corrL1 = Corr(N2+1:end,N2+1:end);
corrL1L2 = Corr(1:N2,N2+1:end); % between-area correlation

corrL1 = triu(corrL1,1);
corrL1 = reshape(corrL1,[N1*N1,1]);
idx = find(corrL1~=0);
corrL1 = corrL1(idx);

corrL2 = triu(corrL2,1);
corrL2 = reshape(corrL2,[N2*N2,1]);
idx = find(corrL2~=0);
corrL2 = corrL2(idx);

corrL1L2 = reshape(corrL1L2,[N1*N2,1]);

% Plotting
edges = -1:0.01:1;

figure()
[bincountsL1,edges] = histcounts(corrL1,edges);
histogram('BinCounts', bincountsL1/sum(bincountsL1), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['Layer 1 - tauIdecay=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL1))])
box off
set(gca,'TickDir','out')

figure()
[bincountsL2,edges] = histcounts(corrL2,edges);
histogram('BinCounts', bincountsL2/sum(bincountsL2), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['Layer 2 - tauIdecay=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL2))])
box off
set(gca,'TickDir','out')

figure()
[bincountsL1L2,edges] = histcounts(corrL1L2,edges);
histogram('BinCounts', bincountsL1L2/sum(bincountsL1L2), 'BinEdges', edges,'FaceColor','black');
xlabel('Spike count correlation')
xlim([-1 1])
ylim([0 0.06])
title(['L1-L2 - tauIdecay=',num2str(paramVal),' in L',num2str(layerOfChange),'; mean=',num2str(mean(corrL1L2))])
box off
set(gca,'TickDir','out')
