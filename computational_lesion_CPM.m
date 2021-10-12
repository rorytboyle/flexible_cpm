% Demo script for applying computational lesions to connectivity matrices
% prior to running CPM
% Author: Rory Boyle
% Date: 12/10/2021

% Take connectivity matrices as used in CPM_Neuromethods_analysis_script.m
% connectivity matrices should be read in as all_mats

%% 1) Apply computational lesion to a functional network
% Download network labels from following url:
% https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv

% Read in network labels (Note: assumes you have downloaded the file to C:\)
network_labels = 'C:\shen_268_parcellation_networklabels.csv'

% Select network to be lesioned
% 1 = Medial frontal network; 2 = Frontoparietal network;
% 3 = Default mode network; 4 = Subcortical-cerebellum network;
% 5 = Motor network; 6 = Visual I network; 7 = Visual II network;
% 8 = Visual association network
network_to_lesion = 5  % lesion motor network

% Apply computational lesion to functional network (e.g. motor network)
lesioned_network_mats = apply_computational_network_lesion_CPM(...
    all_mats, network_to_lesion, network_labels);

% Now use lesioned_network_mats as input in CPM analysis

%% 2) Apply computational lesion to a specific region of interest
clearvars -except all_mats

% Download roi/anatomical labels from following url:
% https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx

% Read in ROI labels (Note: assumes you have downloaded the file to C:\)
ROI_labels = 'C:\1-s2.0-S1053811919309577-mmc2.xlsx'

% Select region to be lesioned, note: region must be listed in the .xlsx
% file containing the labels, otherwise no lesion will be applied
ROI_to_lesion = {'Cerebellum'};  % lesion cerebellum

% Apply computational lesion to ROI (e.g. cerebellum)
lesioned_ROI_mats = apply_computational_ROI_lesion_CPM(all_mats, ROI_to_lesion, ROI_labels);

%% 2) Apply computational lesion to multiple specific regions of interest
clearvars -except all_mats

% Download roi/anatomical labels from following url:
% https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx

% Read in ROI labels (Note: assumes you have downloaded the file to C:\)
ROI_labels = 'C:\1-s2.0-S1053811919309577-mmc2.xlsx'

% Select regions to be lesioned, note: regions must be listed in the .xlsx
% file containing the labels, otherwise no lesion will be applied
ROI_to_lesion = {'Cerebellum', 'Brainstem'};  % lesion cerebellum & brainstem

% Apply computational lesion to ROI (e.g. cerebellum & brainstem)
lesioned_ROI_mats = apply_computational_ROI_lesion_CPM(all_mats, ROI_to_lesion, ROI_labels);
