function [lesioned_mats] = apply_computational_ROI_lesion_CPM(...
    all_mats, ROI_to_lesion, ROI_labels)
% Applies a computational lesion to a connectivity matrix such that nodes 
% from a specified ROI(s) is/are removed from the connectivity matrix.
% The lesioned connectivity matrix can then be used in CPM to derive the 
% importance of the lesioned ROI(s) for prediction of a given phenotype.
% Can also be used to remove nodes from a specific ROI affected by issues
% such as poor scanner coverage. For example, if many scans in a dataset 
% have poor coverage of the cerebellum for a large number of participants, 
% then this can be used to remove all nodes from the cerebellum. 
% ROI labels found here
% https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv
% Nodes are assigned an anatomical label using the Tailarach Atlas as 
% in Salehi et al. (2020) Neuroimage: 
% ROI labels are listed in Table S1 of Salehi et al. (2020), url here:
% https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx

% 
% INPUT:
% all_mats =                (array) m * m * n array where m = number of 
%                           nodes in connectivity matrix and n = number of 
%                           participants. Contains connectivity matrices 
%                           for all participants. Assumes symmetrical 
%                           connectivity matrices.
% ROI_to_lesion =           (cell array) cell array of strings containing 
%                           ROIs to be removed from the connectivity matrix. 
%                           See description above for further info on the 
%                           ROI labels. For example, if you want to remove
%                           all nodes within the cerebellum, then
%                           ROI_to_lesion = {'Cerebellum'}.
% ROI_labels =              (string) file path for .csv file containing
%                           ROI labels: https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx
%
% OUTPUT:
% lesioned_mats =           (array) L * L * n array where L = number of 
%                           nodes in lesioned connectivity matrix and n = 
%                           number of participants. Contains connectivity 
%                           matrices for all participants. 
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 11/10/2021
%
% Read in ROI labels info
ROI_info = importdata(ROI_labels);

% Get anatomical labels
anat_labels = ROI_info.textdata(:,2);

% Get indices of all nodes that are NOT within the ROI(s) to be lesioned
keep_ix = find(~contains(anat_labels(3:end), ROI_to_lesion));

% Remove nodes within lesioned network from connectivity matrices
lesioned_mats = all_mats(keep_ix, keep_ix, :);

end

