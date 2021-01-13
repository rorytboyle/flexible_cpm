function [cpm_target, final_ppts, ix_ppts_to_keep] = ...
    prep_target_CPM(target_csv, file_order, save_path)
% Reads in .csv file containing target variable to predict with 
% connectome-based predictive modeling using YaleMRRC CPM code. Prepares 
% target variable for use in CPM code.
%
% INPUT:
% target_csv =      (string) full file path for directory containing .csv 
%                   files in which each participant's connectivity matrix 
%                   is saved.
% file_order =      n * 1 array containing participant ids for participants
%                   who have data in cpm_predictors (obtained by running
%                   prep_predictors_CPM.m).
% save_path =       (string) file path for output directory. A .mat file
%                   'cpm_predictor_variables', containing cpm_predictors
%                   and file_order will be saved here.
%
% OUTPUT:
% cpm_target =      (double) n * 1 array, containing target variable values,
%                   where n = number of participants. Values will be in the
%                   same order as in final_ppts (i.e. cpm_target(1) will
%                   contain the target variable for the participant in
%                   final_ppts(1)).
% final_ppts =      (double) n * 1 array, containing participant codes,
%                   where n = number of participants. Values will be in the
%                   same order as in cpm_target.
% ix_ppts_to_keep = indices of participants (i.e. connectivity matrices) to
%                   be retained for CPM analysis.
%
% Example usage:    [cpm_target, final_ppts, ix_ppts_to_keep] = ...
%                   prep_predictors_CPM('C:\data\target_variable',...
%                   file_order, 'C:\cpm_data')
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 04/01/2021
%% 1) Prepare participant ids from file_order
% Convert file_order to double
ppt_codes = str2double(file_order);

%% 2) Extract target variable values for participants in final prepped data
% (i.e. for those ppts with connectivity matrices in cpm_predictors)
% Read in .csv file containing target variable
target_all = importdata(target_csv);
target_all = target_all.data;

% extract target variables
[val, pos] = intersect(target_all(:,1), ppt_codes);
target = target_all(pos, :);

final_ppts = target(:, 1);
cpm_target = target(:, 2);

%% 3) Get indices of participants with target variables in ppt_codes
[~, ix_ppts_to_keep] = intersect(ppt_codes, target);

%% 4) Save output
cd(save_path)
save('cpm_target_variable.mat',...
    'cpm_target', 'final_ppts', 'ix_ppts_to_keep')
end