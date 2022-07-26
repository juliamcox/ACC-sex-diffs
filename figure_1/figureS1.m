
% Plots Supplemental Figure 1A, B, C, D

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
basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data'); % where to save data 
cutoff          = inf; % upper cutoff for latency
lowCutoff       = 0; % lower cutoff for latency 
qFile           = 'qLearn_session_all_2022.mat'; %name of file with q-values
fext            = 'fig1';
weightSpreadsheet = fullfile(whereAreWe('figurecode'), 'raw_data','sessionWeights.xlsx');
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort? 
intervalThresh  = 300;
%% Which functions to run



extractFlag(1) = 0; %(1)plot latency distribution by outcome(supp fig 1B,C);
plotFlag(1)    = 0;   

extractFlag(2) = 0; %(2)fit inverse gaussian (supp fig 1A-B) 
plotFlag(2) = 1;


%% Concatenate all data 

if extractAllData
    extractData(aids_m,aids_f, perfThresh,qFile,weightSpreadsheet,ext_laser,binNum,binNum_choice,cohort_opto,sessionLength,intervalThresh);
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
else
    try
        % load all behavior data
        load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
    catch
        extractData(aids_m,aids_f, perfThresh,qFile,weightSpreadsheet,ext_laser,binNum,binNum_choice,cohort_opto,sessionLength,intervalThresh);
        load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
    end
end

behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);
%% Plot


if plotFlag(1) 
    stats_hist=logLatencyHistograms_ctrl(behaviorTable,lowCutoff,cutoff,aids_f,aids_m);
end

if plotFlag(2)
    stats_invGauss= fitLatencyDistributions_inverseGaussian_control(behaviorTable,'none',valType,latencyType);
end

