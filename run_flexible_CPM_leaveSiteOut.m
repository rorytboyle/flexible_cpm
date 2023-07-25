function [behav_pred_pos, behav_pred_neg, behav_pred_combined,...
    parameters_pos, parameters_neg, parameters_combined,...
    pos_mask_all, neg_mask_all, no_node, no_covars] = ...
    run_flexible_CPM_leaveSiteOut(all_behav, all_mats, all_covars, ...
    num_sites, sites, thresh_type, thresh, adjust_stage)
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
% num_sites =       (double) number specifying number of sites in
%                   dataset in order to run leave-site-out cross-validation.
% sites =           (double) n * 1 array, contains the scan sites for each
%                   participant, the first scan site should start from 1.                  
%                    
% thresh_type =     (string) specifies type of feature selection threshold,
%                   must be either 'p-value' or 'sparsity'
% thresh =          (double) If thresh_type = 'p-value, thresh specifies 
%                   p-value for feature selection threshold. If thresh_type
%                   = 'sparsity', thresh specifies % of most highly
%                   correlated edges to be retained (i.e. % threshold for
%                   sparsity threshold).
% adjust_stage      (string) specifies stage of CPM at which covariates
%                   should be adjusted for. 'relate' will adjust for 
%                   covariates at the feature selection step (step 3) by 
%                   running a partial correlation between edges and target 
%                   variable, controlling for covariates. 'fit' will adjust
%                   for covariates at the model fitting step (step 6) by
%                   including covariates in a multiple regression
%                   predicting target variable from network strength.
%                   'both' will adjust for covariates at both feature
%                   selection (step 3) and model fitting (step 6) steps. If
%                   covariates are not being considered, specify an empty 
%                   string ''.
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
% Author: Yihe Weng & Rory Boyle
% Contact: wengy@tcd.ie & rorytboyle@gmail.com
% Date: 24/01/2021
% Last updated: 14/06/2021 added functionality to threshold based on edge
% sparsity (step 3).
% Edited Date:  19/09/2022
%
%% 1) Prepare cross-validated CPM
if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
    % preallocate arrays
    [no_sub, no_node, ~, behav_pred_pos, behav_pred_neg,...
        behav_pred_combined, ~, ~, ~, pos_mask_all, neg_mask_all] =...
        CPM_prep_arrays(all_mats, all_covars, num_sites, adjust_stage);

    parameters_pos = zeros(num_sites, max(sites)+width(all_covars)); % 2 + max(sites) - 2(dummy variable and leave site)
    parameters_neg = zeros(num_sites, max(sites)+width(all_covars));
    parameters_combined = zeros(num_sites, max(sites)+width(all_covars));
else
    [no_sub, no_node, ~, behav_pred_pos, behav_pred_neg,...
    behav_pred_combined, parameters_pos, parameters_neg,...
    parameters_combined, pos_mask_all, neg_mask_all] =...
    CPM_prep_arrays(all_mats, all_covars, num_sites, adjust_stage);
    
end

% specify cross-validation scheme
% kfold_partition = cvpartition(no_sub, 'KFold', k);
site_covars = zeros(length(sites),max(sites));
for s = 1:max(sites)
    site_covars(sites == s,s) = 1;
end


%% 2) Run cross-validated CPM
% loop through each fold and run CPM
for fold = 1:num_sites
    % print message to update
    fprintf('\n Leaving out fold # %6.3f\n',fold);
    
    % make site covars 
    site_lso_covars = site_covars; 
    site_lso_covars(:,fold) = []; % remove test site 
    site_lso_covars(:,end) = []; % remove the last site (dummy variable)
    lso_all_covars = [all_covars,site_lso_covars];
    no_covars = width(lso_all_covars);
    
    % divide data into training and test sets (Step 2 - Shen et al. 2017)
    [ix_train, ix_test, train_behav, train_mats, train_vcts, ...
        train_covars, test_behav, test_mats, test_covars] = ...
        CPM_cv_split_lsocv(all_behav, all_mats, lso_all_covars, no_covars,fold,sites);
     
    % feature selection - relate edges to target variable (Step 3 - Shen et
    % al. 2017)
    % use partial correlation if specified - otherwise use normal
    % correlation
    if strcmp(adjust_stage, 'relate') | strcmp(adjust_stage, 'both')
        [r_mat, p_mat] = CPM_fs_relate_partial(train_vcts, train_behav, ...
            train_covars, no_node);
    else
        [r_mat, p_mat] = CPM_fs_relate(train_vcts, train_behav, no_node);
    end

    % feature selection - select edges (Step 4 - Shen et al. 2017)
    if strcmp(thresh_type, 'p-value');
        [pos_mask, neg_mask] = CPM_fs_select_pvalue(r_mat, p_mat,...
            thresh, no_node);
        
    else strcmp(thresh_type, 'sparsity');
        [pos_mask, neg_mask] = fs_select_sparsity_CPM(r_mat, p_mat,...
            thresh, no_node);  
    end
    
    % calculate network strength in training set (Step 5 - Shen et al.
    % 2017)
    [train_sumpos, train_sumneg, train_sumcombined] = ...
        CPM_network_strength(train_mats, pos_mask, neg_mask, ix_train);
    
    % fit model on training set (Step 6 - Shen et al. 2017)
    % include covars in multiple regression if specified - otherwise use 
    % simple linear regression
    if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
        [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, train_covars, no_covars, ...
        train_sumpos, train_sumneg, train_sumcombined);
    else
        [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, [], 0, ...
        train_sumpos, train_sumneg, train_sumcombined);
    end
    
    % apply model to test set (Step 7 - Shen et al., 2017)
    % account for covars in model application if adjusted for in Step 6
    if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
        [pred_pos, pred_neg, pred_combined] = ...
            CPM_apply_model(test_mats, test_covars, no_covars, ...
            pos_mask, neg_mask, fit_pos, fit_neg, fit_combined);
    else
        [pred_pos, pred_neg, pred_combined] = ...
            CPM_apply_model(test_mats, [], 0, ...
            pos_mask, neg_mask, fit_pos, fit_neg, fit_combined);
    end
    
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