
clear all
close all
%% Parameters
cohort          = 'all';
binNum          = 5;
binNum_choice   = 11;
zscoreFlag      = 1; 
valType         = {'QChosenDiff'};
valType_choice  = {'QDiff'};
latencyType     = {'trialStart'}; 
highCutoff      = inf;
lowCutoff       = 0; 
ageLow          = 1;
ageHigh         = 71;
perfThresh      = .1;
sideThresh      = .4;

basefilename    = whereAreWe('data');


%% Plot q-values (Figure E4b)
 
% Load data 

load(fullfile(basefilename, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'alpha','beta','side','stay');
load(fullfile(basefilename, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'alpha_m','beta_m','side_m','stay_m');

plotQParams(alpha,alpha_m,beta,beta_m,stay,stay_m,side,side_m)

%% Plot choice x value x sex (figure E4c-d)
% load data
load(fullfile(basefilename, sprintf('choiceXvalue_f_bins%d_perfThresh%s_sideThresh%s_human.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_f','bins');
load(fullfile(basefilename, sprintf('choiceXvalue_m_bins%d_perfThresh%s_sideThresh%s_human.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_m');

humanBandit_valueXage_choice(choice_f,choice_m,valType_choice,binNum_choice)
stats_mdl_choice=choiceXValue_linReg_humanBandit_sexDiff(choice_f,choice_m,bins,valType_choice,zscoreFlag);

%% Plot latency x age scatter (figure E4e)
zscoreFlag = 0;
%load data 
load(fullfile(basefilename, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
load(fullfile(basefilename, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');

[fit_m,fit_f] = humanBandit_ageXlatency_scatter(latency_f,latency_m);


%% Plot latency x value x sex (figure E4f-g)
zscoreFlag = 1;
% load data
load(fullfile(basefilename, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
load(fullfile(basefilename, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d_human.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');

stats = plotLatencyXValue_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag);
stats{1}.postHoc_age = humanBandit_valueXage(latency_f,latency_m,valType,latencyType,binNum,zscoreFlag);


