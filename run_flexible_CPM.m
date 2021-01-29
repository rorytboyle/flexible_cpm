function [behav_pred_pos, behav_pred_neg, behav_pred_combined,...
    parameters_pos, parameters_neg, parameters_combined,...
    pos_mask_all, neg_mask_all, no_node, no_covars] = ...
    run_flexible_CPM(all_behav, all_mats, all_covars, k, thresh)
% Runs connectome-based predictive modelling with cross-validation. Enables
% choice of different k-fold cross-validation schemes (can specify LOOCV by
% calling k = number of participants) and allows for covariates to be
% modelled. Returns model predictions for performance evaluation, selected
% edges for visualisation, and model parameters for external validation.
%
% INPUT:
% all_behav =       (array) n * 1 array where n = number of participants.
%                   Contains target variable for all participants
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
% thresh =          (double) p-value for feature selection threshold.
%
% OUTPUT:
% behav_pred_pos =      (array) n * 1 where n = number of participants.
%                       Contains predicted values from
%                       positive network strength model.                       
% behav_pred_neg =      as above for negative network strength model.
% behav_pred_combined = as above for combined network strength model.  
% parameters_pos =      (array) k * m where k = number of folds in CV 
%                       scheme and m = 2 + number of covariates. Contains 
%                       fitted parameters from positive network strength 
%                       model in each CV fold. Rows correspond to folds.
%                       Col 1 = intercept, Col 2 = network strength slope.
%                       If covariates modelled, Col 3 = slope for covar 1,
%                       Col 4 = slope for covar 2, and so on.
% parameters_neg =      as above for negative network strength model.
% parameters_combined = as above for combined network strength model.
% pos_mask_all =        (array) j * j * k where j = number of nodes and k =
%                       number of folds in CV scheme. Each fold contains
%                       binary matrix with 1 indicating that edge was
%                       selected (i.e. positively correlated with the
%                       target variable) in that fold. 
% neg_mask_all =        as above for negatively correlated edges.
% no_node =             (double) number of nodes in connectivity matrix.
%                       Used to extract selected edges.
% no_covars =           (double) number of covariates. If no covariates
%                       included, no_covars = 0. Used to extract average
%                       model parameters.
% 
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 24/01/2021
%
%% 1) Prepare cross-validated CPM
% preallocate arrays
[no_sub, no_node, no_covars, behav_pred_pos, behav_pred_neg,...
    behav_pred_combined, parameters_pos, parameters_neg,...
    parameters_combined, pos_mask_all, neg_mask_all] =...
    CPM_prep_arrays(all_mats, all_covars, k);

% specify cross-validation scheme
kfold_partition = cvpartition(no_sub, 'KFold', k);

%% 2) Run cross-validated CPM
% loop through each fold and run CPM
for fold = 1:k
    % print message to update
    fprintf('\n Leaving out fold # %6.3f\n',fold);
    
    % divide data into training and test sets (Step 2 - Shen et al. 2017)
    [ix_train, ix_test, train_behav, train_mats, train_vcts, ...
        train_covars, test_behav, test_mats, test_covars] = ...
        CPM_cv_split(all_behav, all_mats, all_covars, no_covars, ....
        kfold_partition, fold);
    
    % feature selection - relate edges to target variable (Step 3 - Shen et
    % al. 2017)
    [r_mat, p_mat] = CPM_fs_relate(train_vcts, train_behav, no_node);

    % feature selection - select edges (Step 4 - Shen et al. 2017)
    [pos_mask, neg_mask] = CPM_fs_select(r_mat, p_mat, thresh, no_node);
    
    % calculate network strength in training set (Step 5 - Shen et al.
    % 2017)
    [train_sumpos, train_sumneg, train_sumcombined] = ...
        CPM_network_strength(train_mats, pos_mask, neg_mask, ix_train);
    
    % fit model on training set (Step 6 - Shen et al. 2017)
    [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, train_covars, no_covars, ...
        train_sumpos, train_sumneg, train_sumcombined);
      
    % apply model to test set (Step 7 - Shen et al., 2017)
    [pred_pos, pred_neg, pred_combined] = ...
        CPM_apply_model(test_mats, test_covars, no_covars, ...
        pos_mask, neg_mask, fit_pos, fit_neg, fit_combined);
    
    % store predictions from current fold
    behav_pred_pos(ix_test) = pred_pos;
    behav_pred_neg(ix_test) = pred_neg;
    behav_pred_combined(ix_test) = pred_combined;
    
    % store fitted paramaters from current fold
    parameters_pos(fold, :) = fit_pos';
    parameters_neg(fold, :) = fit_neg';
    parameters_combined(fold, :) = fit_combined'; 
       
    % store edges selected in current fold
    pos_mask_all(:, :, fold) = pos_mask;
    neg_mask_all(:, :, fold) = neg_mask;
    
end
end