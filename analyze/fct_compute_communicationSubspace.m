% Function which computes the communication subspace between two areas,
% based on Semedo et al (2019)

% INPUTS:
% activity1                 : source activity, should be of dimension 
%                             (# observations)x(# variables)
% activity2                 : target activity, should be of dimension
%                             (# observations)x(# variables)
% numDimsUsedForPrediction  : vector containing the numbers of predictive 
%                             dimensions to be tested (eg [0:nNeurons])
% cvNumFolds                : number of cross-validation folds
% bool_print                : 0 (= do not print anything in terminal); or
%                             1 (= print stuff on terminal)
% bool_plot                 : 0 (= do not plot anything); or
%                             1 (= plot figures)

% OUTPUTS:
% optDimReducedRankRegress  : optimal rank for the reduced-rank regression
%                             matrix, as defined by the minimal dimension
%                             for which the prediction performance is
%                             within 1 SEM of the maximal prediction
%                             performance
% optLoss                   : loss when using the optimal rank for the
%                             reduced-rank regression matrix
% optPredPerf_mean          : mean prediction performance over all
%                             cross-validation folds
% optPredPerf_sem           : SEM of prediction performance over all
%                             cross-validation folds


function [optDimReducedRankRegress,optLoss,optPredPerf_mean,optPredPerf_sem] = ...
    fct_compute_communicationSubspace(activity1,activity2,...
    numDimsUsedForPrediction,cvNumFolds,bool_print,bool_plot)

SET_CONSTS_general

% Mean number of spikes per timebins for each neuron separately
meanActivity1 = mean(activity1,1);
meanActivity2 = mean(activity2,1);

% Subtract the mean PSTH per neuron
source_data = double(activity1) - repmat(meanActivity1,size(activity1,1),1);
target_data = double(activity2) - repmat(meanActivity2,size(activity2,1),1);


%% Cross-validated Reduced-Rank Regression

% Initialize default options for cross-validation.
cvOptions = statset('crossval');

% If the MATLAB parallel toolbox is available, uncomment this line to
% enable parallel cross-validation.
% cvOptions.UseParallel = true;

% Regression method to be used.
regressMethod = @ReducedRankRegress;

% Auxiliary function to be used within the cross-validation routine (type
% 'help crossval' for more information). Briefly, it takes as input
% the train and test sets, fits the model to the train set and uses it to
% predict the test set, reporting the model's test performance. Here we
% use NSE (Normalized Squared Error) as the performance metric. MSE (Mean
% Squared Error) is also available.
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
    numDimsUsedForPrediction, 'LossMeasure', 'NSE');

% Cross-validation routine.
% cv1 is of size cvNumFolds x length(numDimsUsedForPrediction)
cvl = crossval(cvFun, target_data, source_data, ...
    'KFold', cvNumFolds, ...
    'Options', cvOptions);

% Stores cross-validation results: mean loss and standard error of the
% mean across folds.
cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];

% To compute the optimal dimensionality for the regression model, call
% ModelSelect:
[optDimReducedRankRegress, optLoss] = ModelSelect(cvLoss, numDimsUsedForPrediction);
if bool_print
    disp(['The optimal dimensionality of the (L1-L2) communication subspace is : ',...
        num2str(optDimReducedRankRegress),' with optLoss :',num2str(optLoss)])
end



%% Prediction performance

% Plot Reduced-Rank Regression cross-validation results
if bool_plot
    figure()
    errorbar(numDimsUsedForPrediction, 1-cvLoss(1,:), cvLoss(2,:),...
        'o--', 'Color', COLOR(target,:),...
        'MarkerFaceColor', COLOR(target,:), 'MarkerSize', 10)
    xlabel('Number of predictive dimensions')
    ylabel('Predictive performance')
end

if optDimReducedRankRegress==0 || optLoss >= 1.0
    optPredPerf_mean = 0;
    optPredPerf_sem = 0;
    if bool_print
        disp('Zero predictive performance since there is no communication subspace...')
    end
else
    optPredPerf_mean = 1-cvLoss(1,numDimsUsedForPrediction==optDimReducedRankRegress);
    optPredPerf_sem = cvLoss(2,numDimsUsedForPrediction==optDimReducedRankRegress);
    if bool_print
        disp(['The predictive performance with ',num2str(optDimReducedRankRegress),...
            ' dimensions is ',num2str(optPredPerf_mean),...
            ' (SEM= ',num2str(optPredPerf_sem),')'])
    end
end

end
