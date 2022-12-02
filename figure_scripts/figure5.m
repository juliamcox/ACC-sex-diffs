%%
clear all
close all

%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 

sigVer     = 'choice';
plotVer    = 'fHist';
sigEvents  = 4; % ipsi choice
plotEvents = {'ipsiLeverPress';'contraLeverPress'};
norm       = 'none';
histLen    = [-2 6; -2 6]; 


%% Figure 5b
choiceHist_maleFemale
%% Figure 5c
% male neurons
savename   = 'choice_m';
plotRegression(plotVer, frameRate, savename,plotEvents,norm,histLen);
% female neurons 
savename   = 'choice_f';
plotRegression(plotVer, frameRate, savename,plotEvents,norm,histLen);
%% Figure 5d and extended figure 9c
[pval, chi2stat, posthoc]=plotPerEncoding_maleFemale(recs_f,recs_m,sigVer,frameRate, rawFlag, zscoreFlag,sigEvents);