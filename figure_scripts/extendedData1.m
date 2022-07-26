
% Plots Extended data Figure 1A, B, C, D

clear all
close all

%% Parameters

% female subjects
aids_f = generateAnimalList('ACC_DMS_nphr_female');
aids_f = cat(1,aids_f,generateAnimalList('ACC_DMS_nphr_yfp_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d1_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d2_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_yfp_female')); 
% male subjects
aids_m = generateAnimalList('ACC_DMS_nphr_male');
aids_m = cat(1,aids_m,generateAnimalList('ACC_DMS_nphr_yfp_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d1_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d2_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_yfp_male')); 

cohort = 'all';

ext             = 'TAB'; % experiment (control sessions)
sessionLength   = 'long';% session length (alternative: 'short';'both')

binNum          = 4; %number of quantile bins for latency
binNum_choice   = 9; %number of quantile bins for choice
zscoreFlag      = 0; %zscore latency?
valType         = 'qChosenDiff'; %which value to plot (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
valType_choice  = 'qDiff';
latencyType     = 'trialInit_thresh'; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press), 'withdrawal' (nose poke entry to exit) 
perfThresh      = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename    = whereAreWe('data'); % where to save data 
cutoff          = inf; % upper cutoff for latency
lowCutoff       = 0; % lower cutoff for latency 
qFile           = 'qLearn_session_all.mat'; %name of file with q-values
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort? 
intervalThresh  = 300;



%% Load data

load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);

%% Figure E1a-b
stats_invGauss= fitLatencyDistributions_inverseGaussian_control(behaviorTable,latencyType);

%% Figure E1c-d
stats_hist=logLatencyHistograms_ctrl(behaviorTable,lowCutoff,cutoff,aids_f,aids_m);

%% Figure E1e-f
latencyXweight(behaviorTable,binNum,aids_f,aids_m)

%% Figure E1g
latencyXtime(behaviorTable,basefilename,binNum,zscoreFlag,perfThresh,qFile,cohort)

%% Figure E1h
load(fullfile(basefilename,['ctrlLatencyStaySwitch_f_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));
load(fullfile(basefilename,['ctrlLatencyStaySwitch_m_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));

stats=plotStaySwitch_ctrl(latency_f,latency_m,{latencyType});
