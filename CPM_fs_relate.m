function [r_mat, p_mat] = CPM_fs_relate(train_vcts, train_behav, no_node)
% Relates functional connectivity in each edge in connectivity matrix to 
% target variable for CPM feature selection.
%
% INPUT: 
% train_vcts =          as described in CPM_cv_split.
% train_behav =         as described in CPM_cv_split.
% no_node =             as described in CPM_prep_arrays.
%
% OUTPUT:
% r_mat =               no_node * no_node array containing Pearson's r
%                       value for correlation of functional
%                       connectivity in each edge with target variable.
% p_mat =               no_node * no_node array containing p-value for
%                       Pearson's r correlation of functional
%                       connectivity in each edge with target variable.
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
% Relate functional connectivity to behaviour
[r_mat,p_mat] = corr(train_vcts',train_behav);

% Reshape from vectors to matrices
r_mat = reshape(r_mat,no_node,no_node);
p_mat = reshape(p_mat,no_node,no_node);
end