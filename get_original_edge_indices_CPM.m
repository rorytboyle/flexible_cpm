function [pos_edges_orig, neg_edges_orig] = ...
    get_original_edge_indices_CPM(pos_edges, neg_edges, node_indices)
% Maps row and column indices from a restricted connectivity matrix (i.e.
% a connectivity matrix with nodes/ROIs removed due to missing coverage) to
% their row and column indices in the full original connectivity matrix.
% This was written specifically for Shen 268 atlas.
%
% INPUT: 
% pos_edges =       (array) p * 4 array where p = number of edges selected 
%                   in positive network model. Col 1 = number of folds in 
%                   which edge selected, Col 2 = linear index of edge in 
%                   vectorised restricted connectivity matrix, Col 3 = row
%                   index of edge in restricted connectivity matrix, Col 4
%                   = col index of edge in connectivity matrix. Returned by
%                   create_masks_CPM.
% neg_edges =       as above for negative network strength model.
% node_indices =    (array) 268*2 array where Col 1 = original node indices
%                   in full connectivity matrix and Col 2 = node indices in
%                   restricted connectivity matrix. Col 2 maps onto Col 1
%                   such that Col 1 contains the original node index (i.e.
%                   in the full connectivity matrix) for the node in the 
%                   same row in Col 2.
%
% OUTPUT:
% pos_edges_orig =  (array) p * 3 array where p = number of edges selected
%                   in positive network model. Col 1 = number of folds in
%                   which edge was selected, Col 2 = row index of edge in
%                   full/original connectivity matrix, Col 3 = col index of 
%                   edge in full/original connectivity matrix.
% neg_edges_orig =  as above for negative network strength.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 26/01/2021
%
%% 1) Account for python zero-indexing - make node_indices 1-indexed
mat_node_indices = node_indices;
% add 1 to node indices (except where new index = 0 as this indicates that 
% this node was not retained in new index i.e. it was a cerebellum or 
% brainstem node)
% cant use simple "node_indices+1" as that adds 1 where node is not in new
% indices (i.e. for cerebellar and brainstem nodes)
for row = 1:length(mat_node_indices)
    % first row will have zero so add 1 here
    if row==1
        mat_node_indices(row,:) = mat_node_indices(row,:)+1;
    end
    if row > 1
        mat_node_indices(row,1) = mat_node_indices(row,1)+1;
        if mat_node_indices(row,2) > 0
            mat_node_indices(row,2) = mat_node_indices(row,2)+1;
        end
    end
end

new_nodes = mat_node_indices(:,2);

%% 2) Get indices in original connectivity matrix for positive network
% get "new/restricted" indices i.e. indices in 205*205 connectivity matrix
current_pos_rows = pos_edges(:,3);
current_pos_cols = pos_edges(:,4);

% map new (i.e. restricted) indices to original indices - positive network
[~, orig_pos_rows] = ismember(current_pos_rows, new_nodes);
[~, orig_pos_cols] = ismember(current_pos_cols, new_nodes);

% add to edge info array
pos_edges_orig = [pos_edges(:,1) orig_pos_rows orig_pos_cols];

%% 3) Get indices in original connectivity matrix for negative network
% get "new/restricted" indices i.e. indices in 205*205 connectivity matrix
current_neg_rows = neg_edges(:,3);
current_neg_cols = neg_edges(:,4);

% map new (i.e. restricted) indices to original indices - positive network
[~, orig_neg_rows] = ismember(current_neg_rows, new_nodes);
[~, orig_neg_cols] = ismember(current_neg_cols, new_nodes);

% add to edge info array
neg_edges_orig = [neg_edges(:,1) orig_neg_rows orig_neg_cols];

end

