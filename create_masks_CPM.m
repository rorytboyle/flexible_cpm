function [pos_edge_mask, neg_edge_mask] = ...
    create_masks_CPM(pos_edges, neg_edges, no_node)
% Creates binary matrices that can be written out as .txt files in 'ascii'
% format for visualisation in BioImageSuite Web. Matrices contain a 1 where
% an edge was selected in positive/negative network strength models.
%
% INPUT: 
% pos_edges =       (array) p * 3 array where p = number of edges selected
%                   in positive network model. Col 1 = number of folds in 
%                   which edge selected, Col 2 = row index of edge in 
%                   connectivity matrix, Col 3 = col index of edge in 
%                   connectivity matrix.
% neg_edges =       as above for negative network strength.
% no_node =         as described in CPM_prep_arrays
%
% OUTPUT:
% pos_edge_mask =   (array) k * k where k = number of nodes in connectivity
%                   matrix. Contains 1 where edge was selected, and 0 where
%                   edge was not selected.
% neg_edge_mask =   as above for negative network strength.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 25/01/2021
%
% create binary matrices for visualisation
pos_edge_mask = zeros(no_node, no_node);
neg_edge_mask = zeros(no_node, no_node);

if ~isempty(pos_edges)  % will ensure mask of zeros returned if no edges selected
    % add 1s to positive binary masks where positive edges were selected
    for ix = 1:size(pos_edges,1)
        pos_edge_mask(pos_edges(ix, 2), pos_edges(ix, 3)) = 1;
    end
end

if ~isempty(neg_edges)  % will ensure mask of zeros returned if no edges selected
% add 1s to negative binary masks where negative edges were selected
    for ix = 1:size(neg_edges,1)
        neg_edge_mask(neg_edges(ix, 2), neg_edges(ix, 3)) = 1;
    end
end
end