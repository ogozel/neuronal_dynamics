% Function that selects a disc of neurons in a grid, with one specific
% neuron in the center, and a given radius.
% It also takes care of the periodic boundary conditions.

% INPUTS:
% centerIdx     : linear index of the neurons that corresponds to the center 
%                 of the disc
% radius        : radius of the disc, in terms of number of neurons
% gridHeight    : height of the grid, in terms of number of neurons
% gridWidth     : width of the grid, in terms of number of neurons
% bool_plot     : if 1: plot the selected disc; if 0: no plot

% OUTPUT:
% idx           : linear indices of the selected neurons


function [idx] = fct_select_discNeurons(centerIdx,radius,gridHeight,gridWidth,bool_plot)

% Grid with all subscripts
[allRows,allCols] = ind2sub([gridHeight,gridWidth],1:gridHeight*gridWidth);

% Subscripts of the center neuron
[rowCenter,colCenter] = ind2sub([gridHeight,gridWidth],centerIdx);

% Euclidian distances between the center neuron and all neurons in the grid
distToCenter = sqrt( (rowCenter-allRows).^2 + (colCenter-allCols).^2 );

% Linear indices of neurons which are within 'radius' of the center neuron
idx = find(distToCenter<=radius);

% Take care of the boundary conditions
boolTB = 0;
boolBB= 0;
boolLB = 0;
boolRB= 0;
if rowCenter-radius < 1 % top boundary
    distToCenterTB = sqrt( (rowCenter-(allRows-gridHeight)).^2 + (colCenter-allCols).^2 );
    idxAddedTB = find(distToCenterTB<=radius);
    idx = cat(2,idx,idxAddedTB);
    boolTB = 1;
end
if rowCenter+radius > gridHeight % bottom boundary
    distToCenterBB = sqrt( ((rowCenter-gridHeight)-allRows).^2 + (colCenter-allCols).^2 );
    idxAddedBB = find(distToCenterBB<=radius);
    idx = cat(2,idx,idxAddedBB);
    boolBB = 1;
end
if colCenter-radius < 1 % left boundary
    distToCenterLB = sqrt( (rowCenter-allRows).^2 + (colCenter-(allCols-gridWidth)).^2 );
    idxAddedLB = find(distToCenterLB<=radius);
    idx = cat(2,idx,idxAddedLB);
    boolLB = 1;
end
if colCenter+radius > gridWidth % right boundary
    distToCenterRB = sqrt( (rowCenter-allRows).^2 + ((colCenter-gridWidth)-allCols).^2 );
    idxAddedRB = find(distToCenterRB<=radius);
    idx = cat(2,idx,idxAddedRB);
    boolRB = 1;
end


% Take care of the corners
distToCenterB = sqrt( (rowCenter-(boolBB*gridHeight)-(allRows-(boolTB*gridHeight))).^2 +...
    (colCenter-(boolRB*gridWidth)-(allCols-boolLB*gridWidth)).^2 );
idxAddedBB = find(distToCenterB<=radius);
idx = unique(cat(2,idx,idxAddedBB));


% Plot the circle of selected neurons
if bool_plot
    tmp = zeros(gridHeight*gridWidth,1);
    tmp(idx) = 1;
    figure()
    imagesc(reshape(tmp,gridHeight,gridWidth))
    xlabel('Neuron index')
    ylabel('Neuron index')
    pbaspect([1 1 1])
    title(['Neuronal disc radius = ',num2str(radius),' neurons'])
end

end