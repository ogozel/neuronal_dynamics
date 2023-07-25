% Function that computes the average subpopulation activity

% INPUTS:
% s0: spike times and indices of neurons which spike
% t_start: beginning of the timewindow we are interested in (ms)
% t_end: end of the timewindow we are interested in (ms)
% idxSelNeur: linear indices of the selected neurons
% bool_plot: 1 to plot the average firing rate per selected neurons, 0
%            otherwise

% OUTPUTS:
% popSpk: average firing rate (Hz) over the whole population
% timeBins: corresponding timebins (ms)


function [popSpk, timeBins] = fct_subpopulationActivity(s0,t_start,t_end,idxSelNeur,bool_plot)

dta = 2; % bin size (ms)
timeBins = t_start+dta/2:dta:t_end-dta/2;
numNeur = numel(idxSelNeur); % number of neurons in selected subpopulation

numframes = numel(timeBins);

popSpk = NaN(numNeur,numframes);

for n=1:numNeur
    
    % Keep only the spikes of the selected neuron
    tmpidx = ismember(s0(2,:),idxSelNeur(n));
    this_s0 = s0(:,tmpidx);
    
    for i=1:numframes
        
        % Find spikes in this time bin
        theseSpk = this_s0(1,:)>timeBins(i)-dta/2 & this_s0(1,:)<=timeBins(i)+dta/2 & this_s0(2,:)>0;
        popSpk(n,i) = sum(theseSpk);
        
    end
end

if bool_plot
    gridHeight = 200;
    gridWidth = 200;
    tmp = zeros(gridHeight*gridWidth,1);
    durationInSec = (t_end-t_start)/1000;
    tmp(idxSelNeur) = sum(popSpk,2)/durationInSec; % average firing rate [Hz]
    figure()
    imagesc(reshape(tmp,gridHeight,gridWidth))
    xlabel('Neuron index')
    ylabel('Neuron index')
    pbaspect([1 1 1])
end


end