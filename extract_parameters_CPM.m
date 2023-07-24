function [int_pos_ntwrk, int_neg_ntwrk, int_combined_ntwrk,...
    slope_pos_ntwrk, slope_neg_ntwrk, slope_combined_ntwrk,...
    slope_pos_covars, slope_neg_covars, slope_combined_covars] = ...
    extract_parameters_CPM(parameters_pos, parameters_neg,...
    parameters_combined, adjust_stage, no_covars)
% Calculates average model parameters (i.e. regression coefficients) across
% folds that can then be applied to independent datasets for external
% validation.
%
% INPUT: 
% parameters_pos =      as described in run_flexible_CPM.
% parameters_neg =      as described in run_flexible_CPM.
% parameters_combined = as described in run_flexible_CPM.
% adjust_stage =        as described in run_flexible_CPM.
% no_covars =           as described in CPM_prep_arrays. 
%
% OUTPUT:
% int_pos_ntwrk =           (double) intercept for positive network 
%                           strength model.
% int_neg_ntwrk =           as above for negative network strength.
% int_combined_ntwrk =      as above for combined network strength.
% slope_pos_ntwrk =         (double) slope for positive network strength.
% slope_neg_ntwrk =         as above for negative network strength.
% slope_combined_ntwrk =    as above for combined network strength.
% slope_pos_covars =        (array) 1 * p where p = number of covariates.
%                           Contains slope for covariates in positive
%                           network strength model. Col 1 = average slope
%                           for covariate 1, col 2 = average slope for
%                           covariate 2, and so on.
% slope_neg_covars =        as above for negative network strength.
% slope_combined_covars =   as above for combined network strength.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 25/01/2021
% Updated: 01/06/2021
% Updated: 14/07/2023
%
% Get average intercepts for each network strength model across folds
int_pos_ntwrk = mean(parameters_pos(:,1));
int_neg_ntwrk = mean(parameters_neg(:,1));
int_combined_ntwrk = mean(parameters_combined(:,1));

% Get average slopes for network strength in each model across folds
slope_pos_ntwrk = mean(parameters_pos(:,2));
slope_neg_ntwrk = mean(parameters_neg(:,2));
slope_combined_ntwrk = mean(parameters_combined(:,2));

% Get average slopes for covariates in each model across folds
if no_covars > 0 && (strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both'))
    slope_pos_covars =  mean(parameters_pos(:,3:end));
    slope_neg_covars = mean(parameters_neg(:,3:end));
    slope_combined_covars = mean(parameters_combined(:,3:end));
    
% this block necessary for general analysis script to run correctly as
% covar parameters will require some value (e.g. NaN)
elseif no_covars > 0 && strcmp(adjust_stage, 'relate')
    slope_pos_covars = nan(1,no_covars);
    slope_neg_covars = nan(1,no_covars);
    slope_combined_covars = nan(1,no_covars);
    
else
    slope_pos_covars = [];
    slope_neg_covars = [];
    slope_combined_covars = [];
end
end