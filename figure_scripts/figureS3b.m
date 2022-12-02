clear all
close all

%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 
qFile      = 'qLearn_session_all_2022.mat'; %name of file with q-values

% parameters for loading behavior file
binNum = 4; % number of value bins
perfThresh = 0.1; % performance threshold 
sessionLength = 'both'; % which sessions to include
binNum_choice = 9; % number of value bins for choice
intervalThresh = 300; % interval threshold for session 
cohort_opto = 'ACC_DMS_nphr';

[r,p]=imaging_control_indVariability(perfThresh,qFile,sessionLength,binNum,binNum_choice,intervalThresh,cohort_opto,recs_f,recs_m);