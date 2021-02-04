function [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, train_covars, no_covars, ...
        train_sumpos, train_sumneg, train_sumcombined)
% Fits linear regression model where network strength values (and
% covariates if specified) are related to the target variable.
%
% INPUT:
% train_behav =         as described in CPM_cv_split.
% train_covars =        as described in CPM_cv_split.
% no_covars =           as described in CPM_prep_arrays.
% train_sumpos =        as described in CPM_network_strength.
% train_sumneg =        as described in CPM_network_strength.
% train_sumcombined =   as described in CPM_network_strength.
%
% OUTPUT:
% fit_pos =             1 * p array where p = 2 + no_covars. Stores
%                       parameters from fitted model for positive
%                       network. Row 1 contains intercept, row contains
%                       slope for network strength. If covars included, row
%                       3 contains slope for covar 1, row 4 contains slope 
%                       for covar 2, and so on.
% fit_neg =             same as above for negative network.
% fit_combined =        same as above for combined network.
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
% Create independent variables for each network strength model - add column
% of ones to network strength values
x_pos = [ones(size(train_behav)) train_sumpos];
x_neg = [ones(size(train_behav)) train_sumneg];
x_combined = [ones(size(train_behav)) train_sumcombined];

% Add covariates to independent variables for each network strength model
if no_covars>0
    for covar = 1:no_covars
        x_pos(:,covar+2) = train_covars(:,covar);
        x_neg(:,covar+2) = train_covars(:,covar);
        x_combined(:,covar+2) = train_covars(:,covar);
    end
end

% Fit each network strength model
fit_pos = regress(train_behav, x_pos);
fit_neg = regress(train_behav, x_neg);
fit_combined = regress(train_behav, x_combined);
end