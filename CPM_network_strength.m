function [train_sumpos, train_sumneg, train_sumcombined] = ...
        CPM_network_strength(train_mats, pos_mask, neg_mask, ix_train)
% Calculates network strength values for positive, negative, and combined
% networks. Network strength = summed connectivity of edges that are
% related to target variable. Combined network strength is calculated using
% the method described by Greene et al. 2018 Nature Communications
% https://doi.org/10.1038/s41467-018-04920-3 where combined network strength
% is the difference between positive and negative network strength.
%
% INPUT:
% train_mats =          as described in CPM_cv_split.
% pos_mask =            as described in CPM_fs_select.
% neg_mask =            as described in CPM_fs_select.
% ix_train =            as described in CPM_cv_split.
%
% OUTPUT:
% train_sumpos =        n * 1 array where n = number of participants in
%                       current training set. Contains summed positive
%                       network strength values (i.e. summed connectivity
%                       of edges that are positively related to target
%                       variable).
% train_sumneg =        n * 1 array where n = number of participants in
%                       current training set. Contains summed negative
%                       network strength values (i.e. summed connectivity
%                       of edges that are negatively related to target
%                       variable).
% train_sumcombined =   n * 1 array where n = number of participants in
%                       current training set. Contains difference in
%                       network strength between positive and negative 
%                       networks. Defined according to Greene et al. 2018 
%                       Nature Communications as positive network minus
%                       negative network.
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
% Create arrays to store network strength values
train_sumpos = zeros(length(ix_train),1);
train_sumneg = zeros(length(ix_train),1);

% Get sum of all positive and negative thresholded edges in training set
% participants. Divide by 2 control for the fact that matrices are
% symmetric.
for ss = 1:size(train_sumpos)
    train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
    train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
end

% Calculate combined network strength
train_sumcombined = train_sumpos - train_sumneg;
end
