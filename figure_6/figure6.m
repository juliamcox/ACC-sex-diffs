%%Figure 6 
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
sideVer    = 'v26'; % version of regression for relative side value
chosenVer  = 'v10'; % version of regression for relative chosen value
sigLevel   = 0.01;
%% 6a 
% load summary of all neurons
load(fullfile(whereAreWe('bilinear'),dataset,chosenVer,'summary.mat'),'pvals','pvals_m','thisFnames','thisFnames_m');

plotExamples_bilinearReg(pvals,sigLevel,chosenVer,thisFnames);
plotExamples_bilinearReg(pvals_m,sigLevel,chosenVer,thisFnames_m);

%% 6b
bilinearRegression_summary(chosenVer,dataset,rawFlag)

%% 6c
outcomeIdx = 5;
bilinearRegression_valueSummary(chosenVer,dataset,linspace(1,100,4),outcomeIdx,rawFlag,frameRate,qFile)
%% 6d plot individual neurons' fluorescence by relative side value

% load summary of all neurons
load(fullfile(whereAreWe('bilinear'),dataset,sideVer,'summary.mat'),'pvals','pvals_m','thisFnames','thisFnames_m');

plotExamples_bilinearReg(pvals,sigLevel,sideVer,thisFnames);
plotExamples_bilinearReg(pvals_m,sigLevel,sideVer,thisFnames_m);

%% 6e
bilinearRegression_summary(sideVer,dataset,rawFlag)

%% 6f
leverIdx = 3;
bilinearRegression_valueSummary(sideVer,dataset,linspace(0,100,4),leverIdx,rawFlag,frameRate,qFile)