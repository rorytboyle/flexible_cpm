function [cpm_predictors, file_order] = prep_predictors_CPM(folder,...
    n, save_path)
% Reads in .csv files containing z-scored connectivity matrices in order
% to run connectome-based predictive modeling using YaleMRRC CPM code.
%
% INPUT:
% folder =          (string) file path for directory containing .csv files
%                   in which each participant's connectivity matrix is
%                   saved.
% n =               (double) number of participants with connectivity
%                   matrices. Must match number of .csv files in 'folder'.
% save_path =       (string) file path for output directory. A .mat file
%                   'cpm_predictor_variables', containing cpm_predictors
%                   and file_order will be saved here.
%
% OUTPUT:
% cpm_predictors =  (double) x * x * n array, containing x * x connectivity
%                   matrices for all n participants, where x = number of 
%                   ROIs and n = number of participants.
% file_order =      (cell) n * 1 array containing file names. Position of
%                   file name in file_order corresponds to position in 
%                   cpm_predictors
%
% Example usage:    [cpm_predictors] = prep_data_CPM(...
%                   'C:\connectivity_matrices', 350', 'C:\cpm_data')
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 04/01/2021
%
%% 1) Ensure correct number of files specified
% read in files
cd(folder);
files = dir('*.csv');

% Check number of .csv files == n
if n ~= length(files)
   error("'n' (number of participants) is not equal to number of .csv files in 'folder'")
end

%% 2) Load connectivity matrices
% create cell for storing file names
file_order = cell(n, 1);

% create array for storing all connectivity matrices
sample_conn_mx = load(files(1).name);
cpm_predictors = zeros(length(sample_conn_mx), length(sample_conn_mx), n);

% loop through files and load connectivity matrices
for i = 1:n
    file_name = files(i).name;
    conn_mx = load(file_name);
    
    file_order{i} = file_name;
    cpm_predictors(:, :, i) = conn_mx;
end

%% 3) Save output
cd(save_path)
save('cpm_predictor_variables_all_ppts.mat', 'cpm_predictors', 'file_order')
end