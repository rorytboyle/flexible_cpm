% This analysis script can be used to run computational lesions to 
% connectivity matrices of a specific network or multiple networkts to see 
% how networks contribute to the CPM prediction 
% Authors: Rory Boyle & Yihe Weng
% Emails: rorytboyle@gmail.com & wengy@tcd.ie
% Date: 07/10/2021
% Date edited: 20/05/2022
clear; clc;
%% 1) load data and prepare variables
% load predictor variables (i.e. connectivity matrices)
load('W:\AOMIC\output\kfold_10\cpm_predictor_variables_all_ppts.mat','cpm_predictors')
all_mats = cpm_predictors;

% load target data
load('W:\AOMIC\output\kfold_10\cpm_target_variable.mat','cpm_target');
all_behav = cpm_target;

% load covariates
covars_path = 'W:\AOMIC\doc\covariates.csv';
data = readtable(covars_path);
covar_names = {'age','sex'};
all_covars = table2array(data(:, covar_names));

% Specify output file paths
output_path = 'W:\AOMIC\output\kfold_10\lesionCPM';

% load oringal CPM parameters
load('W:\AOMIC\output\kfold_10\IQ_CPM.mat','model_info');

%% 2) Apply computational lesion to a functional network
% Note: Assumes the Shen atlas regions/labels is downloaded to C:\
% The path of Shen atlas regions/labels can be changed in the function:
% extract_lesion_network_mats_CPM.m
% Download network labels from following url:
% https://www.nitrc.org/frs/download.php/8072/shen_268_parcellation_networklabels.csv

% Download roi/anatomical labels from following url:
% https://ars.els-cdn.com/content/image/1-s2.0-S1053811919308249-mmc1.xlsx  

% Method:'1' means using network index to remove the network's connectivity 
% if method is 1, network_to_lesion should be the index of lesion network
% for example: network_to_lesion = 5; 
% method: '0' means using roi/anatomical labels to remove the network's 
% connectivity (could be a specific region or multiple regions). for example: 
% network_to_lesion = {'Cerebellum'}; or network_to_lesion = {'Cerebellum', 'Brainstem'};
method = '1';
% Select network to be lesioned
% 1 = Medial frontal network; 2 = Frontoparietal network;
% 3 = Default mode network; 4 = Subcortical-cerebellum network;
% 5 = Motor network; 6 = Visual I network; 7 = Visual II network;
% 8 = Visual association network
%     network_to_lesion = 5; 
network_to_lesion = 5; % lesion motor network
% network_to_lesion = {'Cerebellum'};  % lesion cerebellum
lesioned_network_mats = extract_lesion_network_mats_CPM(all_mats,method,network_to_lesion);

%% 3) Specify model inputs
adjust_stage = model_info.adjust_stage;
k = model_info.k;
iterations = model_info.iterations;  
thresh_type = model_info.fs_thresh_type;
thresh = model_info.fs_thresh;
freq_thresh = model_info.freq_thresh;
model_info.permutation_test = 'no'; %run_permutation = 'no';
model_name = model_info.model_name;
model_info.output_path = output_path;

%% 4) Preallocate arrays for storing results and parameters 
% Preallocate arrays for storing CPM predicted values and results
[behav_pred_pos_all, behav_pred_neg_all, behav_pred_combined_all, ...
    R_pos, P_pos, R_neg, P_neg, R_combined, P_combined,...
    rsq_pos, rsq_neg, rsq_combined, mae_pos, mae_neg, mae_combined] = ...
    prep_results_arrays_CPM(lesioned_network_mats, iterations);

% Preallocate arrays for storing CPM parameters
if strcmp(adjust_stage, 'relate')
    [~, ~, int_pos_ntwrk, int_neg_ntwrk,...
    int_combined_ntwrk, slope_pos_ntwrk, slope_neg_ntwrk, ...
    slope_combined_ntwrk, ~, ~, ~] = prep_parameters_arrays_CPM(...
    lesioned_network_mats,all_covars, k, iterations);   
else 
    [~, ~, int_pos_ntwrk, int_neg_ntwrk,...
    int_combined_ntwrk, slope_pos_ntwrk, slope_neg_ntwrk, ...
    slope_combined_ntwrk, slope_pos_covars, slope_neg_covars,...
    slope_combined_covars] = prep_parameters_arrays_CPM(...
    lesioned_network_mats,all_covars, k, iterations);   
end

%% 5) Run CPM, evaluate model performance, extract selected edges + model parameters
parfor i = 1:iterations
    fprintf('\n Running iteration # %6.3f\n',i);
     
    % Run CPM
     [behav_pred_pos_all(:,i), behav_pred_neg_all(:,i), ...
         behav_pred_combined_all(:,i), parameters_pos, parameters_neg,...
         parameters_combined, pos_mask, neg_mask, no_node, no_covars] = ...
         run_flexible_CPM(all_behav, lesioned_network_mats, all_covars, k, ...
         thresh_type, thresh, adjust_stage); 
    
    % Evaluate model performance
    [R_pos(i), P_pos(i), R_neg(i), P_neg(i), R_combined(i), P_combined(i),...
        rsq_pos(i), rsq_neg(i), rsq_combined(i), mae_pos(i), mae_neg(i),...
        mae_combined(i)] = evaluate_CPM(all_behav, behav_pred_pos_all(:,i),...
        behav_pred_neg_all(:,i), behav_pred_combined_all(:,i));

    % Extract model parameters
    if strcmp(adjust_stage, 'relate')
        [int_pos_ntwrk(i), int_neg_ntwrk(i), int_combined_ntwrk(i),...
            slope_pos_ntwrk(i), slope_neg_ntwrk(i), slope_combined_ntwrk(i),...
            ~, ~, ~] = extract_parameters_CPM(parameters_pos, parameters_neg,...
            parameters_combined, adjust_stage, no_covars); 
    else
        [int_pos_ntwrk(i), int_neg_ntwrk(i), int_combined_ntwrk(i),...
            slope_pos_ntwrk(i), slope_neg_ntwrk(i), slope_combined_ntwrk(i),...
            slope_pos_covars(i,:), slope_neg_covars(i,:), ...
            slope_combined_covars(i,:)] = ...
            extract_parameters_CPM(parameters_pos, parameters_neg,...
            parameters_combined,adjust_stage, no_covars); 
    end 

% Store positive and negative edges masks - accounting for parallel loop
    c_pos_mask_all{i,1} = pos_mask;
    c_neg_mask_all{i,1} = neg_mask;  
end


% Add positive and negative edge masks
pos_mask_all = [];
neg_mask_all = [];
for j = 1:iterations
    pos_mask_all = cat(3, c_pos_mask_all{j,1}, pos_mask_all);
    neg_mask_all = cat(3, c_neg_mask_all{j,1}, neg_mask_all);
end

%% 6) Extract selected edges 
% set number of folds across all iterations of model
k_all = k*iterations;

% no_node not available after parfor loop so call again here
no_node = size(lesioned_network_mats,1);

% get indices of edges in each network that are correlated with the
% target variable in >= number of total folds specified by freq_thresh
[pos_edges_all, neg_edges_all, pos_edges_thresh, neg_edges_thresh] = ...
    extract_edges_CPM(pos_mask_all, neg_mask_all, no_node, k_all, freq_thresh);

%% 7) Store parameters, predicted values, edges, and results in structs and save data
if strcmp(adjust_stage, 'relate')
    parameters = struct('int_pos_ntwrk', int_pos_ntwrk, 'int_neg_ntwrk', ...
    int_neg_ntwrk, 'int_combined_ntwrk', int_combined_ntwrk,...
    'slope_pos_ntwrk', slope_pos_ntwrk, 'slope_neg_ntwrk', slope_neg_ntwrk,...
    'slope_combined_ntwrk', slope_combined_ntwrk,...
    'pos_mask_all', pos_mask_all, 'neg_mask_all', neg_mask_all); 
else
    parameters = struct('int_pos_ntwrk', int_pos_ntwrk, 'int_neg_ntwrk', ...
    int_neg_ntwrk, 'int_combined_ntwrk', int_combined_ntwrk,...
    'slope_pos_ntwrk', slope_pos_ntwrk, 'slope_neg_ntwrk', slope_neg_ntwrk,...
    'slope_combined_ntwrk', slope_combined_ntwrk, 'slope_pos_covars',...
    slope_pos_covars, 'slope_neg_covars', slope_neg_covars,...
    'slope_combined_covars', slope_combined_covars,...
    'pos_mask_all', pos_mask_all, 'neg_mask_all', neg_mask_all);
end

predictions = struct('pos_preds_all_folds', behav_pred_pos_all, ...
    'neg_preds_all_folds', behav_pred_neg_all, 'combined_preds_all_folds',...
    behav_pred_combined_all);

edges = struct('pos_edges_all', pos_edges_all, 'neg_edges_all', neg_edges_all,...
    'pos_edges_thresh', pos_edges_thresh, 'neg_edges_thresh', neg_edges_thresh);

results = struct('R_pos', R_pos, 'R_neg', R_neg, 'R_combined', R_combined,...
    'rsq_pos', rsq_pos, 'rsq_neg', rsq_neg, 'rsq_combined', rsq_combined,...
    'mae_pos', mae_pos, 'mae_neg', mae_neg', 'mae_combined', mae_combined);

% Get mean parameters, predicted values, and results across iterations
if iterations > 1
    mean_parameters = structfun(@mean, parameters, 'UniformOutput', false);
    mean_predictions = structfun(@(x) mean(x, 2), predictions, 'UniformOutput', false);
    mean_results = structfun(@mean, results, 'UniformOutput', false);
end

% save data
variable_list = {'model_info', 'parameters', 'predictions', 'results', 'edges', ....
    'mean_parameters', 'mean_predictions', 'mean_results'};

save_file = [model_info.output_path filesep model_name '.mat'];
for variable_ix = 1:length(variable_list)
    if exist(variable_list{variable_ix})
        if exist(save_file)
            save(save_file, variable_list{variable_ix}, '-append');
        else
            save(save_file, variable_list{variable_ix});
        end
    end
end
