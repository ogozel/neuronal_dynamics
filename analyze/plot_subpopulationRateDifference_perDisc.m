% Analyze the results from 'compute_subpopulationRate.m'
% (1) Plot the average firing rate over all neurons as a function of time
%       for each realization separately
% (2) Plot the integral of the square of the difference in average firing
%       rate for discs as a function of disc radius

% Run this scripts on the local machine after having saved the results from
% running 'compute_subpopulationRate.m' on the HPC cluster locally.

clear all
%close all
clc

dataFolder = '../data_cluster_tauIdecay/';
% '../data_cluster_sigmaI/'
% '../data_cluster_tauIdecay/'
filename1 = 'popRate_seed1000_V1_tau24_s2_1to21s';
filename2 = 'popRate_seed1000_V2_tau24_s2_1to21s';
% 'popRate_seed1000_V1_sigma0d1_s2_1to21s'
%       'popRate_seed1000_V2_sigma0d1_s2_1to21s'
% 'popRate_seed1000_V1_sigma0d3_s2_1to21s'
%       'popRate_seed1000_V2_sigma0d3_s2_1to21s'
% 'popRate_seed1000_V1_tau24_s2_1to21s'
%       'popRate_seed1000_V2_tau24_s2_1to21s'


bool_smooth = 0; % smooth the firing rate using a boxcar filter or not
bool_only50 = 0; % randomly select only 50 neurons per disc


% Timewindow for the spike counts
binSize = 2; % [ms]

% Total number of excitatory neurons
NE = 40000;


%% Load data
data1 = load([dataFolder, filename1, '.mat']);
data2 = load([dataFolder, filename2, '.mat']);

% Center of spike count bins (same for both realizations)
binCenters = data1.timeBins;

%% Smooth firing rate

if bool_smooth
    
    WindowLength = 100; %6*20;  % 200 ms
    % GaussianSmoothSigma = 20; % 100 ms
    % GaussianSmoothAlpha = (WindowLength-1)/(2*GaussianSmoothSigma);
    % wind = gausswin(WindowLength,GaussianSmoothAlpha);
    wind = boxcar(WindowLength);
    smWin = wind./sum(wind); % normalize the convolution kernel to 1
    
    spkMat = data1.popSpk';
    spkMat = reshape(spkMat, 1, size(spkMat,1), size(spkMat,2));
    Master_smoothedBa = fct_smoothPSTH(spkMat,binCenters,smWin,WindowLength,binSize);
    spk1 = squeeze(Master_smoothedBa)';
    
    spkMat = data2.popSpk';
    spkMat = reshape(spkMat, 1, size(spkMat,1), size(spkMat,2));
    Master_smoothedBa = fct_smoothPSTH(spkMat,binCenters,smWin,WindowLength,binSize);
    spk2 = squeeze(Master_smoothedBa)';
    
else
    
    spk1 = data1.popSpk;
    spk2 = data2.popSpk;
    
end

%% Plot the difference in average firing rate 
% as a function of time between the two realizations

% data.popSpk is a (# neurons) x (# timebins) array, trials of 20s, so
% 10,000 timebins of 2ms (and 40,000 excitatory neurons)

%%%  over all neurons %%%
avgFR1 = mean(spk1, 1) * (1000/binSize);
avgFR2 = mean(spk2, 1) * (1000/binSize);

figure()
subplot(2,1,1)
plot(binCenters, avgFR1)
hold on
plot(binCenters, avgFR2)
xlim([binCenters(1) binCenters(end)])
ylabel('Pop act [spk/s]')
box off

subplot(2,1,2)
plot(binCenters, avgFR1-avgFR2)
xlabel('Timebin (ms)')
ylabel('Delta Pop act [spk/s]')
xlim([binCenters(1) binCenters(end)])
ylim([-45 45])
box off

%%%  over the neurons sampled from a small disc %%%
% Select a disc of neurons
radius = 5;
centerIdx = randperm(NE); % Center of the disc
idx = fct_select_discNeurons(centerIdx, radius, sqrt(NE), sqrt(NE), 0);

avgFR1 = mean(spk1(idx,:), 1) * (1000/binSize);
avgFR2 = mean(spk2(idx,:), 1) * (1000/binSize);

figure()
subplot(2,1,1)
plot(binCenters, avgFR1)
hold on
plot(binCenters, avgFR2)
xlim([binCenters(1) binCenters(end)])
ylabel('Pop act [spk/s]')
box off

subplot(2,1,2)
plot(binCenters, avgFR1-avgFR2)
xlabel('Timebin (ms)')
ylabel('Delta Pop act [spk/s]')
xlim([binCenters(1) binCenters(end)])
box off
set(gca, 'TickDir', 'out')
ylim([-60 60])


%% Plot the integral of the square of the difference in average firing
%       rate for discs as a function of disc radius

% Total number of draws
nDraws = 10;

Master_radius = [5:19, 20:20:100]; %[1:10, 20:20:100];
Master_results = NaN(length(Master_radius),nDraws);

for r=1:length(Master_radius)
    
    % Radius of disc of neurons
    radius = Master_radius(r);
    
    for d=1:nDraws
        
        % Center of the disc
        centerIdx = randperm(NE);
        
        % Select a disc of neurons
        idx = fct_select_discNeurons(centerIdx, radius, sqrt(NE), sqrt(NE), 0);
        
        if bool_only50
            % Pick only 50 neurons out of this disc
            idx = idx(randperm(length(idx), 50));
        end
        
        % Average firing rate in the disc of neurons in each realization
        % [Hz]
        avgFR1_disc = mean(spk1(idx,:), 1) * (1000/binSize);
        avgFR2_disc = mean(spk2(idx,:), 1) * (1000/binSize);
        
        %%% sqrt(N) * mean[(avgFR1-avgFR2)^2] / mean[avgFR^2] %%%
        
        % Absolute value of the difference [Hz]
        sq_diff = (avgFR1_disc - avgFR2_disc).^2;
        
        % Mean (over the two realizations of 20s) average (over the 
        % neurons) firing rate
        avgFR = mean(cat(2, avgFR1_disc, avgFR2_disc));
                
        % Save the results
        Master_results(r,d) = sqrt(length(idx)) * mean(sq_diff) / avgFR^2;
        
    end
end

% Mean
whole_mean = mean(Master_results,2);
whole_std = std(Master_results,1,2); %/sqrt(nDraws);

% Plot the mean square of the difference
figure(1)
hold on
errorbar(Master_radius, whole_mean, whole_std)
xlabel('Disc radius')
%ylabel('Mean square of the difference')

figure(2)
hold on
semilogy(Master_radius, whole_mean)
xlabel('Disc radius')
set(gca, 'YScale', 'log')

figure(13)
hold on
x = Master_radius;
y1 = transpose(whole_mean - whole_std);
y2 = transpose(whole_mean + whole_std);
X = [x, fliplr(x)];        % create continuous x value array for plotting
Y = [y1, fliplr(y2)];      % create y values for out and then back
fill(X,Y,'y','FaceAlpha',0.3); 
set(gca,'TickDir','out')


