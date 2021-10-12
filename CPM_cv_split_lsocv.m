function [ix_train, ix_test, train_behav, train_mats, train_vcts, ...
    train_covars, test_behav, test_mats, test_covars] = ...
    CPM_cv_split_lsocv(all_behav, all_mats, all_covars, no_covars,fold,s_index)
% This cv split code is for parallelised CPM with leave site out cross-validation.

% INPUT: 
% no_covars =           number of covariates. If no covariates included,
%                       no_covars = 0
% kfold_partition =     cvpartition object that defines random
%                       nonstratified partition of data for k-fold CV
% fold =                fold number in current loop
% other args =          as described in run_flexible_CPM.
% s_index =             scan sites index
%
% OUTPUT:
% ix_train =            indices of training set participants
% ix_test =             indices of test set participant(s)
% train_behav =         target variable for training set
% train_mats =          connectivity matrices for training set
% train_vcts =          p * no_subs array containing vectors for each
%                       participant. Each vector contains connectivity data
%                       reshaped from matrix to vector (to enable
%                       correlation with target variable for feature
%                       selection)
% train_covars =        covariates for training set
% test_behav =          target variable for test set
% test_mats =           connectivity matrices for test set  
% test_covars =         covariates for test set
%
% Author: Yihe Weng & Rory Boyle
% Contact: wengy@tcdie & rorytboyle@gmail.com
% Date: 12/10/2021
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
%% 1) Get indices of training and test set participants
ix_train = find(s_index ~= fold);
ix_test = find(s_index == fold);

%% 2) Remove test set participants from training set data
train_mats = all_mats(:, :, ix_train);
train_vcts = reshape(train_mats,[],size(train_mats,3));

train_behav = all_behav(ix_train,:);

% leave out subject from covars
if no_covars>0
    train_covars = all_covars(ix_train,:);
else
    train_covars = [];
end

%% 3) Get data for test set participants
test_mats = all_mats(:, :, ix_test);
test_behav = all_behav(ix_test,:);
if no_covars>0
    test_covars = all_covars(ix_test,:);
else
    test_covars = [];
end
end