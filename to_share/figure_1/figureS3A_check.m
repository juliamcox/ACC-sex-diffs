% Plots Supplemental Figure 2A

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
valType         = {'qTot'}; %which value to plot (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
valType_choice  = {'qDiff'};
latencyType     = {'trialStart_thresh'}; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press), 'withdrawal' (nose poke entry to exit) 
perfThresh      = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data'); % where to save data 
cutoff          = inf; % upper cutoff for latency
lowCutoff       = 0; % lower cutoff for latency 
qFile           = 'qLearn_session_all.mat'; %name of file with q-values
bandwidth       = .5; % bandwidth for kernel density function
ptiles          = [20:20:100]; % 
fext            = 'fig1';
%% Which functions to run

extractFlag(1) = 0; % [(1)extract latencies and choice in value bins and by previous outcome (Figure 1C,H,B,G); 
plotFlag(1)    = 1; 



%% Extract data

if extractFlag(1)
     ctrlXvalueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext);
end



%% Plot latency x value x sex 

if plotFlag(1)
    % try to load extracted value, if not saved, extract
    try
        load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    catch
        ctrlXvalueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,'fig1');        
        load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    end
    stats_latency = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins);
end

