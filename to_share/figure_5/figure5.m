%%Figure 5 
clear all
close all

%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 
qFile      = 'qLearn_session_all_2022.mat'; %name of file with q-values

%% Figure 5a

%% Figure 5b
choiceScatter_maleFemale(recs_f,recs_m,frameRate,0,0)
%% Figure 5c
sigVer     = 'choice';
plotVer    = 'fHist';
savename   = 'choice';
sigEvents  = 4; % ipsi choice
plotEvents = {'ipsiLeverPress';'contraLeverPress'};
norm       = 'none';
histLen    = [-2 6; -2 6]; 
% male neurons 
plotRegression(recs_m, sigVer, plotVer, frameRate, rawFlag, zscoreFlag, savename, sigEvents,plotEvents,norm,histLen,0.01);
% female neurons 
plotRegression(recs_f, sigVer, plotVer, frameRate, rawFlag, zscoreFlag, savename, sigEvents,plotEvents,norm,histLen,0.01);
%% Figure 5d and figure S13c
ver = 'choice';
sigEvents = 4;
[pval, chi2stat, posthoc]=plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,sigEvents);