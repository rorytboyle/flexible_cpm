function [behav_pred_pos_all, behav_pred_neg_all, behav_pred_combined_all,...
    R_pos, P_pos, R_neg, P_neg, R_combined, P_combined,...
    rsq_pos, rsq_neg, rsq_combined, mae_pos, mae_neg, mae_combined] = ...
    prep_results_arrays_CPM(all_mats, iterations)
% Prepares arrays for storing results of multiple iterations of CPM with
% k-fold cross-validation or single iteration of CPM with LOOCV.
%
% INPUT:
% all_mats =                (array) m * m * n array where m = number of 
%                           nodes in connectivity matrix and n = number of
%                           participants. Contains connectivity matrices 
%                           for all participants. Assumes symmetrical
%                           connectivity matrices.
% iterations =              (double) number of iterations of CPM with 
%                           k-fold CV.
%
% OUTPUT:
% behav_pred_pos_all =      (array) n*iterations array of zeros.
% behav_pred_neg_all =      (array) n*iterations array of zeros.
% behav_pred_combined_all = (array) n*iterations array of zeros.
% R_pos =                   (array) iterations*1 array of zeros.
% P_pos =                   same as above.
% R_neg =                   same as above.
% P_neg =                   same as above.
% R_combined =              same as above.
% P_combined =              same as above.
% rsq_pos =                 same as above.
% rsq_neg =                 same as above.
% rsq_combined =            same as above.
% mae_pos =                 same as above.
% mae_neg =                 same as above.
% mae_combined =            same as above.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 29/04/2021
% Updated: 31/05/2021 added arrays to store predictions from each
% iteration of k-fold CV

% Get number of participants
n = size(all_mats, 3);

% Get number of nodes
n_nodes = size(all_mats, 1);

% Prep arrays for storing predicted values
behav_pred_pos_all = zeros(n, iterations);
behav_pred_neg_all = zeros(n, iterations);
behav_pred_combined_all = zeros(n, iterations);

% Prep arrays for storing results
R_pos = zeros(iterations,1);
P_pos = zeros(iterations,1);
R_neg = zeros(iterations,1);
P_neg = zeros(iterations,1);
R_combined = zeros(iterations,1);
P_combined = zeros(iterations,1);
rsq_pos = zeros(iterations,1);
rsq_neg = zeros(iterations,1);
rsq_combined = zeros(iterations,1);
mae_pos = zeros(iterations, 1);
mae_neg = zeros(iterations, 1);
mae_combined = zeros(iterations,1);

end

