clear all
%% Parameters

recs_m     = generateAnimalList('ACC_DMS_imaging_male');
recs_f     = generateAnimalList('ACC_DMS_imaging_female');

frameRate  = 20;
rawFlag    = 1; 
zscoreFlag = 1; 

%% Figure 4b (plots fluorescence traces for each neuron 1 by 1 with the spatial profile from CNMFe for a randomly selected time) 
plotExampleTraces(recs_m{5})

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
%% Figure 4f and Extended Data 9b
ver = 'outcome';
sigEvents = 6;
[pval, chi2stat, posthoc]=plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,sigEvents);
%% Figure 4h and Extended Data 9d
ver = 'stay';
sigEvents = [5,7];
[pval_stay, chi2stat_stay, posthoc_stay]=plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,sigEvents);

%% Extended data figure 9a 
histLen    = [-2 6; 0 8; -2 6; -2 6; 0 8; 0 8]; 
plotVer    = 'fHist';
norm       = 'none';

plotEvents = {'nosePokeEntry';'leverPresentation';'LLeverPress';'RLeverPress';'rewardDelivery';'CSMinus'};
plotRegression(plotVer,frameRate, ['Basic3'],plotEvents,norm,histLen);