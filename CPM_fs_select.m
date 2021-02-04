function [pos_mask, neg_mask] = CPM_fs_select(r_mat, p_mat, thresh, no_node)
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
% Date: 24/01/2021
%
% NOTE: This code is a modification of codeshare_behavioralprediction.m
% which was written by Xilin Shen and Emily Finn (as obtained here: 
% https://www.nitrc.org/frs/download.php/8071/nn_code.zip).
% Copyright 2015 Xilin Shen and Emily Finn 
% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 
%
% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.
% 
% Create arrays to store edges
pos_mask = zeros(no_node,no_node);
neg_mask = zeros(no_node,no_node);

% Find edges with positive correlation to target variable that survive
% threshold
pos_edges = find(r_mat > 0 & p_mat < thresh);
% Find edges with positive correlation to target variable that survive
% threshold
neg_edges = find(r_mat < 0 & p_mat < thresh);

% Store edges in masks
pos_mask(pos_edges) = 1;
neg_mask(neg_edges) = 1;
end