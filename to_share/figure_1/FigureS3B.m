% Figure 2 D-F, Supp Fig 5, Supp Fig 6 

clear 
close all

%% Parameters
cohort          = 'ACC_DMS_nphr';
ext             = 'LAS'; 
 
binNum          = 4; % number of bins for plotting latency
binNum_choice   = 9; % number of bins for plotting choice
zscoreFlag      = 0; % zscore latencies?
epochs          = {'PrevOutcome'}; % which laser epochs to plot
%epochs_extract  = {'PrevOutcome';'PrevNP';'ITI';'NP';'PrevITI'};
valType         = {'qTot'}; % value type for latency
valType_choice  = {'qDiff'}; % value type for choice
latencyType     = {'trialInit_thresh'}; % trialStart (trial initiation latency), leverPress (lever press latency), withdraw (nose poke exit latency)
perfThresh      = .1; % performance threshold (difference between reward and no reward stay probability)
qFile           = 'qLearn_session_all_2022.mat'; % file with value estimates from Q-learning model
% latency cutoffs
cutoff          = inf; 
lowCutoff       = 0; 
 

basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data');

fext            = 'figS3';

load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');



%% Plot mean latencies in value quantile bins
% Figure 1E, SuppFig 5B
laserType      = 1;
stats_latencyQuant = latencyQuant_opto(cohort,behaviorTable,valType{nv},binNum,laserType,zscoreFlag,latencyType{1});


