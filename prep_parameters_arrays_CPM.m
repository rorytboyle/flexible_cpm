function [pos_mask_all_iterations, neg_mask_all_iterations,...
    int_pos_ntwrk_all, int_neg_ntwrk_all, int_combined_ntwrk_all,...
    slope_pos_ntwrk_all, slope_neg_ntwrk_all, slope_combined_ntwrk_all, ...
    slope_pos_covars_all, slope_neg_covars_all, slope_combined_covars_all]...
    = prep_parameters_arrays_CPM(all_mats, all_covars, k, iterations)
% Prepares arrays for storing parameters obtained from multiple iterations 
% of CPM with k-fold cross-validation
%
% INPUT:
% all_mats =        (array) m * m * n array where m = number of nodes in
%                   connectivity matrix and n = number of participants.
%                   Contains connectivity matrices for all participants.
%                   Assumes symmetrical connectivity matrices.
% all_covars =      (array) n * p array where n = number of participants
%                   and p = number of covariates. If no covariates are to
%                   be entered in model, specify all_covars = [].
% k =               (double) number specifying number of folds in
%                   cross-validation scheme. k = 10 will run 10-fold
%                   cross-validation, k = 5 will run 5-fold, and k = number
%                   of participants will run leave-one-out
%                   cross-validation.
% iterations =      (double) number of iterations of CPM with k-fold CV.
%
% OUTPUT:
% pos_mask_all_iterations =     (array) iterations*1 array of zeros
% neg_mask_all_iterations =     same as above.
% int_pos_ntwrk_all =           (array) iterations*1 array of zeros
% int_neg_ntwrk_all =           same as above.
% int_combined_ntwrk_all =      same as above.
% slope_pos_ntwrk_all =         same as above.
% slope_neg_ntwrk_all =         same as above.
% slope_combined_ntwrk_all =    same as above.
% slope_pos_covars_all =        (array) iterations*(number of covars) array
%                               of zeros. If no covars are specified, an 
%                               empty array will be returned.
% slope_neg_covars_all =        same as above.
% slope_combined_covars_all =   same as above.

% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 28/05/2021
% Updated: 31/08/2021 changed slope_pos_covars_all(and same arrays for neg
% and combined networks) to use zeros(iterations, n_covars) instead of
% zeros(iterations, size(all_covars, 2)) for clarity. Updated description
% of function output.

% Get number of nodes
n_nodes = size(all_mats, 1);

% Get number of covars
n_covars = size(all_covars, 2);

pos_mask_all_iterations = zeros(n_nodes, n_nodes, (k*iterations));
neg_mask_all_iterations = zeros(n_nodes, n_nodes, (k*iterations));
int_pos_ntwrk_all = zeros(iterations,1);
int_neg_ntwrk_all = zeros(iterations,1);
int_combined_ntwrk_all = zeros(iterations,1);
slope_pos_ntwrk_all = zeros(iterations,1);
slope_neg_ntwrk_all = zeros(iterations,1);
slope_combined_ntwrk_all = zeros(iterations,1);
slope_pos_covars_all = zeros(iterations,n_covars);
slope_neg_covars_all = zeros(iterations,n_covars);
slope_combined_covars_all = zeros(iterations,n_covars);

end
