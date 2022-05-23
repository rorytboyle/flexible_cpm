%% Download the functional connectivity for AOMIC dataset and prepare the csv
% file for target, covariable and matrix data for CPM analysis  
% Download the AOMIC data from https://github.com/eminSerin/NBSPredict_SpringerNature
% git clone --recursive git@github.com:eminSerin/NBSPredict_SpringerNature.git
clear;
clc;
%% speficy path and load target and covariables data
% assume the data is downloaded here: 'W:\AOMIC'
rootdir = 'W:\AOMIC\NBSPredict_SpringerNature\data\conn_mat'; 
load('W:\AOMIC\NBSPredict_SpringerNature\data\design_IQ.mat');
N = 20; % subject's number
design = design(1:N,:);
%% make doc,connectivity matrices and output dirtory for CPM analysis
doc_dir = fullfile('W:\AOMIC\doc');   % for saving target,covarible or sites csv files
matrix_dir =  fullfile('W:\AOMIC\connectivity_matrices'); % for saving connectivity matrices
output = fullfile('W:\AOMIC\output'); % output folder
if ~exist(doc_dir,'dir')
    mkdir(doc_dir)
end
if ~exist(matrix_dir,'dir')
    mkdir(matrix_dir)
end
if ~exist(output,'dir')
    mkdir(output)
end     

%% save connectivity csv for N subjects
subfolder = dir(fullfile(rootdir,'sub*'));  
for j = 1: N
    sub_name = subfolder(j).name;
    sub(j) = str2double(sub_name(5:end-4));    
    load(fullfile(rootdir,sub_name));      
    writematrix(data,fullfile(matrix_dir,[sub_name(5:end-4),'.csv']));
    clear data
end

%% save target csv and covariables csv
% make target csv 
T_behav = table();
T_behav.subid = sub';
T_behav.all_behav = design(:,2);
% make covariables csv (age, sex)
T_covar = table();
T_covar.subid = sub';
T_covar.age = design(:,3);
T_covar.sex = design(:,4);
% save csv file
covar_filename = fullfile(doc_dir,'covariates.csv');
behav_filename = fullfile(doc_dir,'behav_data.csv');
writetable(T_covar,covar_filename);
writetable(T_behav,behav_filename);

%% make false site information for leave site out
% all subject id
all_sub = [];
for j = 1: length(subfolder)
    sub_name = subfolder(j).name;
    all_sub(j) = str2double(sub_name(5:end-4));    
end
sites = repmat([1:1:7],1,123);
sites = sites(randperm(length(sites)));
T_sites = table();
T_sites.subid = all_sub';
T_sites.sites = sites';  
T_sites = T_sites(1:N,:);
% save sites csv file
site_filename = fullfile(doc_dir,'sites.csv'); 
writetable(T_sites,site_filename);
