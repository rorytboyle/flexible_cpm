function [pos_edges_all, neg_edges_all, pos_edges_thresh, neg_edges_thresh] = ...
    extract_edges_CPM(pos_mask_all, neg_mask_all, no_node, k_all, freq_thresh)
% Extracts edges that were selected in all folds and in a specified %
% of folds in both positive and network strength models.
%
% INPUT: 
% pos_mask_all =        as described in run_flexible_CPM.
% neg_mask_all =        as described in run_flexible_CPM.
% no_node =             as described in CPM_prep_arrays. 
% k_all =               number of folds across all iterations. If k = n 
%                       (LOOCV) and 1 iteration run, then k_all = n. If k =
%                       10 and 100 iterations run, then k_all = 1000.
% freq_thresh =         (double) frequency threshold for edge selection. If
%                       freq_thresh = 100, edge must be selected (i.e.
%                       significantly correlated with target variable) in
%                       100% (i.e. all) folds of the model.
%
% OUTPUT:
% pos_edges_all =       (array) p * 4 array where p = number of edges
%                       selected in at least 1 fold of positive network 
%                       model. Col 1 = number of folds in which edge
%                       selected, Col 2 = linear index of edge in
%                       vectorised connectivity matrix, Col 3 = row index
%                       of edge in connectivity matrix, Col 4 = col index
%                       of edge in connectivity matrix.
% neg_edges_all =       as above for negative network strength.
% pos_edges_thresh =    (array) p * 4 array where p = number of edges
%                       selected in at least x folds of positive network 
%                       model. X = number of folds determined by 
%                       freq_thresh. Col 1 = number of folds in which edge
%                       selected, Col 2 = linear index of edge in
%                       vectorised connectivity matrix, Col 3 = row index
%                       of edge in connectivity matrix, Col 4 = col index
%                       of edge in connectivity matrix.
% neg_edges_thresh =    as above for negative network strength.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 25/01/2021
% Updated: 14/06/2021 (changed freq_thresh to be decimal percentage i.e.
% between 0 and 1, so that it was coherent with thresh for feature
% selection in fs_select_sparsity_CPM.m and run_flexible_CPM.m
%
% Preallocate arrays
pos_edge_freq = zeros(no_node, no_node);
neg_edge_freq = zeros(no_node, no_node);

% Loop through masks and sum edges across folds (i.e. if k/number of folds
% = 10 and an edge was selected in all folds, then edge will have total
% value of 10).
% RB 31/05/2021 - COULD SPEED THIS BY USING sum(pos_mask_all, 3);
% sum(neg_mask_all,3);
for i = 1:no_node
    for j = 1:no_node
        pos_edge_freq(i,j) = sum(pos_mask_all(i, j, :));
        neg_edge_freq(i,j) = sum(neg_mask_all(i, j, :));
    end
end

% get linear indexes and frequencies for edges selected in at least 1 fold
pos_lin_ix = find(pos_edge_freq>=1);
freqs_pos = pos_edge_freq(pos_lin_ix);
neg_lin_ix = find(neg_edge_freq>=1);
freqs_neg = neg_edge_freq(neg_lin_ix);

% get row and column indices in connectivity matrices for edges selected in
% at least 1 fold
[pos_row, pos_col] = find(pos_edge_freq>=1);
[neg_row, neg_col] = find(neg_edge_freq>=1);

% add to arrays
pos_edges_all = [freqs_pos pos_lin_ix pos_row pos_col];
neg_edges_all = [freqs_neg neg_lin_ix neg_row neg_col];

% calculate required number of folds from frequency threshold
folds_thresh = k_all * freq_thresh;

% get thresholded edges
pos_edges_thresh = pos_edges_all((pos_edges_all(:,1) >= folds_thresh),:);
neg_edges_thresh = neg_edges_all((neg_edges_all(:,1) >= folds_thresh),:);
end
