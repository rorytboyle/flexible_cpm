function [lesioned_mats] = apply_computational_network_lesion_CPM(...
    all_mats, network_to_lesion, network_labels)
% Applies a computational lesion to a connectivity matrix such that all
% nodes from a functional network are removed from the connectivity matrix.
% The lesioned connectivity matrix can then be used in CPM to derive the 
% importance of the lesioned functional network for prediction of a given
% phenotype.
% Network labels found here
% https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv
% Function requires this .csv file
% Network numbers are outlined in Finn et al. (2015) Nature Neuroscience
% https://doi.org/10.1038/nn.4135 and correspond to following functional 
% networks
% #1 = Medial frontal network
% #2 = Frontoparietal network
% #3 = Default mode network
% #4 = Subcortical-cerebellum network
% #5 = Motor network
% #6 = Visual I network
% #7 = Visual II network
% #8 = Visual association network
% 
% INPUT:
% all_mats =                (array) m * m * n array where m = number of 
%                           nodes in connectivity matrix and n = number of 
%                           participants. Contains connectivity matrices 
%                           for all participants. Assumes symmetrical 
%                           connectivity matrices.
% network_to_lesion =       (double) label of network to be removed from 
%                           connectivity matrix. See description above for
%                           further info on the functional network labels.
% network_labels =          (string) file path for .csv file containing
%                           network labels: https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv
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
% Read in network labels
labels = importdata(network_labels).data;

% Get indices of all nodes that are NOT within the network to be lesioned
keep_ix = find(labels(:,2)~=network_to_lesion);

% Remove nodes within lesioned network from connectivity matrices
lesioned_mats = all_mats(keep_ix, keep_ix, :);

end

