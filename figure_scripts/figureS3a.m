% Figure S3a

clear 
close all

%% Parameters
cohort          = 'ACC_DMS_nphr';
ext             = 'LAS'; 


 
binNum          = 4; % number of bins for plotting latency
binNum_choice   = 9; % number of bins for plotting choice
zscoreFlag      = 0; % zscore latencies?

valType         = {'qChosenDiff'}; % value type for latency
valType_choice  = {'qDiff'}; % value type for choice
latencyType     = {'trialInit_thresh'}; % trialInit (trial initiation latency), leverPress (lever press latency)
perfThresh      = .1; % performance threshold (difference between reward and no reward stay probability)
qFile           = 'qLearn_session_all.mat'; % file with value estimates from Q-learning model
% latency cutoffs
cutoff          = inf; 
lowCutoff       = 0; 
sessionLength   = 'long';
intervalThresh  = 300; 
groupingVar     = 'none'; % for fitting inverse gaussian distribution
laserType       = 1;
basefilename    = whereAreWe('data');




%% Load behavior data

load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');


%% Plot control session value modulation v. effect of inhibtion

[corr_coeff,p_corr]=opto_control_indVariability(behaviorTable,cohort,ext,zscoreFlag,groupingVar,laserType);