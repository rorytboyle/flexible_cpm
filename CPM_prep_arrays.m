function [no_sub, no_node, no_covars, behav_pred_pos, behav_pred_neg,...
    behav_pred_combined, parameters_pos, parameters_neg,...
    parameters_combined, pos_mask_all, neg_mask_all] = ...
    CPM_prep_arrays(all_mats, all_covars, k, adjust_stage)
% Prepares and preallocates arrays for connectome-based predictive model.
%
% INPUT: as described in run_flexible_CPM.
%
% OUTPUT:
% no_sub =              number of participants in dataset
% no_node =             number of nodes in connectivity matrix
% no_covars =           number of covariates
% behav_pred_pos =      no_sub * 1 array to store predictions for positive
%                       network
% behav_pred_neg =      as above for negative network.
% behav_pred_combined = as above for combined network.
% parameters_pos =      no_sub * p array where p = 2 + no_covars. Stores
%                       parameters from fitted models using positive
%                       network. Row 1 contains intercept, row contains
%                       slope for network strength. If covars included, row
%                       3 contains slope for covar 1, row 4 contains slope 
%                       for covar 2, and so on.
% parameters_neg =      as above for negative network.
% parameters_combined = as above for combined network.
% pos_mask_all =        no_node * no_node * no_sub array which stores the
%                       edges selected in each fold for the positive
%                       network.
% neg_mask_all =        as above for negative network.
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
% get number of ppts, nodes, arrays
no_sub = size(all_mats,3);
no_node = size(all_mats,1);
no_covars = size(all_covars, 2);

% preallocate arrays to store predictions
behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred_combined = zeros(no_sub, 1);

% preallocate arrays to store paramaters of regression models in each fold
% col 1 = intercept, col 2 = network strength slope.
% if covars are included, and are adjusted for at the model fitting step
% (step 6), then col 3 = covar 1 intercept, col 4 = covar 1 intercept, 
% and so on...
if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
    parameters_pos = zeros(k, 2 + no_covars);
    parameters_neg = zeros(k, 2 + no_covars);
    parameters_combined = zeros(k, 2 + no_covars);
else
    parameters_pos = zeros(k, 2);
    parameters_neg = zeros(k, 2);
    parameters_combined = zeros(k, 2);
end

% preallocate arrays to store selected edges in each fold
pos_mask_all = zeros(no_node, no_node, k);
neg_mask_all = zeros(no_node, no_node, k);

end