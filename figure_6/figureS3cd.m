%%Figure S3 
clear all
close all

%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 
qFile      = 'qLearn_session_all'; %name of file with q-values
dataset    = 'dataset_single1';
totVer    = 'v12'; % version of regression for relative side value
sigLevel   = 0.01;


%% S3c
bilinearRegression_summary(totVer,dataset,rawFlag)

%% S3d
outcomeIdx = 5;
bilinearRegression_valueSummary(totVer,dataset,linspace(1,100,4),outcomeIdx,rawFlag,frameRate,qFile)
