% Function which computes Factor Analysis (FA)

% INPUTS:
% data      : data on which to perform Factor Analysis, 
%             size: (number of neurons)x(number of timepoints)
% zDimList  : list of dimensions to try out for FA
% numFolds  : number of cross-validation folds for FA

% OUTPUTS:
% dim       : results for each dimension from zDimList
% optDimPE  : optimal dimensionality as determined by the minimal
%             prediction error
% optDimLL  : optimal dimensionality as determined by the maximal
%             log-likelihood


function [dim, optDimPE, optDimLL] = fct_computeFA(data,zDimList,numFolds)


% Subtract the average activity per neuron
meanX = mean(data,2);
X_meanSub = data - repmat(meanX,1,size(data,2));

% crossvalidate_fa_OG.m calls fastfa.m, which uses FA by default, but could
% also use PPCA (typ - 'fa' (default) or 'ppca')
[dim, optDimPE, optDimLL] = ...
    crossvalidate_fa_OG(X_meanSub,'zDimList',zDimList,'numFolds',numFolds);


end