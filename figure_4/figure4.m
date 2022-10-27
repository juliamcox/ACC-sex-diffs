clear all
%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 
qFile      = 'qLearn_session_all_2022.mat'; %name of file with q-values

%% Figure 4b (plots fluorescence traces for each neuron 1 by 1 with the spatial profile from CNMFe for a randomly selected time) 
plotExampleTraces(recs_m{5})
%% Figure 4c

%% Figure 4d
outcomeHist_maleFemale
%% Figure 4e
sigVer     = 'outcome';
plotVer    = 'fHist';
sigEvents  = 6; % CS no reward  
plotEvents = {'rewardDelivery';'CSMinus'};
norm       = 'none';
histLen    = [0 8; 0 8]; 
sigLevel   = 0.01;
% male neurons 
plotRegression(plotVer,frameRate, [sigVer '_m'],plotEvents,norm,histLen);
% female neurons 
plotRegression(plotVer,frameRate, [sigVer '_f'],plotEvents,norm,histLen);
%% Figure 4f and S13b
ver = 'outcome';
sigEvents = 6;
[pval, chi2stat, posthoc]=plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,sigEvents);
%% Figure 4g
events = {'outcome'};
sigOutcome_staySwitch_maleFemale(rawFlag,frameRate,zscoreFlag,'Basic3',6,qFile,events,[0 8],'Basic3');
%% Figure 4h and S13d
ver = 'stay';
sigEvents = [5,7];
[pval_stay, chi2stat_stay, posthoc_stay]=plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,sigEvents);
%% Figure S9b
% parameters for loading behavior file
binNum = 4; % number of value bins
perfThresh = 0.1; % performance threshold 
sessionLength = 'both'; % which sessions to include
binNum_choice = 9; % number of value bins for choice
intervalThresh = 300; % interval threshold for session 
cohort_opto = 'ACC_DMS_nphr';
[r,p]=imaging_control_indVariability(perfThresh,qFile,sessionLength,binNum,binNum_choice,intervalThresh,cohort_opto,recs_f,recs_m);
%% Figure S13a
histLen    = [-2 6; 0 8; -2 6; -2 6; 0 8; 0 8]; 
plotEvents = {'nosePokeEntry';'leverPresentation';'LLeverPress';'RLeverPress';'rewardDelivery';'CSMinus'};
plotRegression(cat(1,recs_f,recs_m), 'Basic3', 'fHist', frameRate, rawFlag, zscoreFlag, 'Basic3', 1:6,plotEvents,'zscore',histLen,1);