% Compute the within-area shared dimensionality, by performing Factor
% Analysis to get the shared covariance matrix, and then assessing
% dimensionality (1) through the participation ratio measure, and (2) by
% determining how many dimensions are needed to explain >95% of the shared
% variability.

% INPUTS:
% data              : data on which to compute the shared covariance matrix
%                     size: (number of neurons)x(number of timepoints)
% zDimList          : list of dimensions to try out for FA
% numFolds          : number of cross-validation folds for FA
% bool_forceLow     : if 1, force a low-dimensional representation of 
%                     arbitrary chosen 'numDimForced' dimensions
% numDimForced      : number of dimensions to force the low-dimensional 
%                     representation to
% bool_print        : if 1, display dimensionality on command line                  

% OUTPUTS:
% d_shared_PR       : dimensionality based on Participation Ratio measure
% d_shared_95       : number of dimensions needed to explain >95% of the
%                     variance
% eigenvectors      : eigenvectors of the shared covariance matrix
% eigenvalues       : eigenvalues of the shared covariance matrix
% sharedCovarianceMatrix : shared covariance matrix obtained through Factor
%                          Analysis, size: (num_selNeur x num_selNeur)
% privateCovarianceElements : elements of the diagonal private covariance
%                             matrix obained through Factor Analysis
% L                 : factor loadings, size: (num_selNeur x M)


function [d_shared_PR, d_shared_95, eigenvectors, eigenvalues,...
          sharedCovarianceMatrix, privateCovarianceElements, L] = ...
          fct_compute_sharedDim(data,zDimList,numFolds,bool_forceLow,...
                                numDimForced,bool_print)


% Compute Factor Analysis : X = LL^T + Ph = WW^T + Psi
[dim, ~, optDimLL] = fct_computeFA(data,zDimList,numFolds);

% Optimal dimensionality that maximizes the cross-validated log likelihood
M = optDimLL;

% Elements of the (diagonal) private component of the covariance matrix
privateCovarianceElements = dim(zDimList==M).estParams.Ph;

% Estimate the shared covariance matrix using the optimal dimensionality 
% that maximizes the cross-validated log likelihood (as in Semedo2019)
L = dim(zDimList==M).estParams.L;
if bool_forceLow
    L_low = L(:,1:numDimForced);
    sharedCovarianceMatrix = L_low*L_low';
else
    sharedCovarianceMatrix = L*L';
end

% Eigenvalue decomposition
[eigenvectors,D,~] = svd(sharedCovarianceMatrix);
eigenvalues = diag(D);
eigenvalues = sort(eigenvalues,'descend');

% Choose the direction of the eigenvectors so that most elements are
% positive
for i=1:size(eigenvectors,2)
    tmp = eigenvectors(:,i);
    if (length(find(tmp<0))/length(find(tmp>0))) > 1
        eigenvectors(:,i) = -eigenvectors(:,i);
    end
end

% Percentage of variance in the shared covariance matrix explained
perc_varExpl = cumsum(eigenvalues)/sum(eigenvalues);

% Smallest number of dimensions that capture 95% of the variance in the 
% shared covariance matrix
d_shared_95 = find(perc_varExpl>0.95,1);

% Participation Ratio to estimate within-area dimensionality
d_shared_PR = (sum(eigenvalues)^2)/sum(eigenvalues.^2);

if bool_print
    disp(['Shared dimensionality estimated by Participation Ratio : ',...
        num2str(d_shared_PR),' (out of ',...
        num2str(size(sharedCovarianceMatrix,1)),' neurons)'])
    disp(['Shared dimensionality estimated by number of dim that explain >95% of variance : ',...
        num2str(d_shared_95),' (out of ',num2str(size(sharedCovarianceMatrix,1)),' neurons)'])
end

end