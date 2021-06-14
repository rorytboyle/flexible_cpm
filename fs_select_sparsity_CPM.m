function [pos_mask, neg_mask] = fs_select_sparsity_CPM(r_mat, p_mat,...
    thresh, no_node)
% Selects edges where p-value for correlation between functional 
% connectivity and target variable is below specified threshold.
%
% INPUT: 
% train_vcts =          as described in CPM_fs_relate
% train_behav =         as described in CPM_fs_relate.
% thresh =              as describedin run_flexible_CPM.
% no_node =             as described in CPM_prep_arrays.
%
% OUTPUT:
% pos_mask =            no_node * no_node array where 1 indicates edge was
%                       positively related to target variable in current
%                       fold and correlation survived threshold.
% neg_mask =            no_node * no_node array where 1 indicates edge was
%                       negatively related to target variable in current
%                       fold and correlation survived threshold.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 14/06/2021
%
% Create arrays to store edges
pos_mask = zeros(no_node,no_node);
neg_mask = zeros(no_node,no_node);

% Calculate number of edges to be returned in each network
% Get number of unique mx elements by squaring number of nodes, then
% subtracting number of nodes (i.e. minus diagonal), then divide by 2
% (account for symmetry). Multiply num of unique mx elements by sparsity
% threshold to get the % of edges to be returned.
max_edges = round((((no_node^2) - no_node)/2)*thresh);

% Reshape r_mat to vector
r_vec = reshape(r_mat, [], 1);

% Get min correlation value that survives sparsity threshold (for pos
% network)
min_corr = min(maxk(r_vec, max_edges));

% Get max correlation value that survives sparsity threshold (for neg
% network)
max_corr = max(mink(r_vec, max_edges));

% Find edges with positive correlation to target variable that survive
% threshold
pos_edges = find(r_mat >= min_corr);

% Find edges with negative correlation to target variable that survive
% threshold
neg_edges = find(r_mat <= max_corr);

% Store edges in masks
pos_mask(pos_edges) = 1;
neg_mask(neg_edges) = 1;
end