function lesioned_network_mats = extract_lesion_network_mats_CPM(all_mats,method,network_to_lesion)
% Demo script for applying computational lesions to connectivity matrices
% prior to running CPM
% Authors: Rory Boyle & Yihe Weng
% Emails: rorytboyle@gmail.com & wengy@tcd.ie
% Date: 12/10/2021
% Date edited: 12/05/2022

% all_mat: predictor variables (i.e. connectivity matrices) across all subjects, could be from
% the whole-brain CPM output (cpm_predictor_variables_all_ppts.mat,'cpm_predictors')

% method: 1 -- using network index to remove the lesion-network connectivity
% if method is 1, network_to_lesion should be the index of lesion network
% for example: network_to_lesion = 5;  % lesion motor network

% 1 = Medial frontal network; 2 = Frontoparietal network;
% 3 = Default mode network; 4 = Subcortical-cerebellum network;
% 5 = Motor network; 6 = Visual I network; 7 = Visual II network;
% 8 = Visual association network

% method: 0 -- using network label name to remove the lesion-network
% connectivity (can be a specific region or multiple regions), network_to_lesion 
% should be the ROI name of lesion network. for example: 
% network_to_lesion = {'Cerebellum'}; or network_to_lesion = {'Cerebellum', 'Brainstem'};

% Take connectivity matrices as used in CPM_Neuromethods_analysis_script.m
% connectivity matrices should be read in as all_mats
if strcmp(method,'1')   
    %% 1) Apply computational lesion to a functional network
    % Download network labels from following url:
    % https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv

    % Read in network labels (Note: assumes you have downloaded the file to C:\)
    network_labels = 'C:\shen_268_parcellation_networklabels.csv';

% Select network to be lesioned
% network_to_lesion = 5;  % lesion motor network

    % Apply computational lesion to functional network (e.g. motor network)
    lesioned_network_mats = apply_computational_network_lesion_CPM(...
        all_mats, network_to_lesion, network_labels);

% Now use lesioned_network_mats as input in CPM analysis
else 
%% 2) Apply computational lesion to a specific region of interest

    % Download roi/anatomical labels from following url:
    % https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx

    % Read in ROI labels (Note: assumes you have downloaded the file to C:\)
    ROI_labels = 'C:\1-s2.0-S1053811919308249-mmc1.xlsx';

    % Select region to be lesioned, note: region must be listed in the .xlsx
    % file containing the labels, otherwise no lesion will be applied
%     network_to_lesion = {'Cerebellum'};  % lesion cerebellum

    % Apply computational lesion to ROI (e.g. cerebellum)
    lesioned_network_mats = apply_computational_ROI_lesion_CPM(all_mats, network_to_lesion, ROI_labels);
end 
% %% 2) Apply computational lesion to multiple specific regions of interest
% clearvars -except all_mats
% 
% % Download roi/anatomical labels from following url:
% % https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx
% 
% % Read in ROI labels (Note: assumes you have downloaded the file to C:\)
% ROI_labels = 'C:\1-s2.0-S1053811919309577-mmc2.xlsx'
% 
% % Select regions to be lesioned, note: regions must be listed in the .xlsx
% % file containing the labels, otherwise no lesion will be applied
% ROI_to_lesion = {'Cerebellum', 'Brainstem'};  % lesion cerebellum & brainstem
% 
% % Apply computational lesion to ROI (e.g. cerebellum & brainstem)
% lesioned_ROI_mats = apply_computational_ROI_lesion_CPM(all_mats, ROI_to_lesion, ROI_labels);
