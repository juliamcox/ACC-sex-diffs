% Plots Figures 1B,C,E,G,H


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

ext               = 'TAB'; % experiment (control sessions)
sessionLength     = 'long';% session length (alternative: 'short';'both')

binNum            = 4; %number of quantile bins for latency
binNum_choice     = 9; %number of quantile bins for choice
zscoreFlag        = 0; %zscore latency?
valType           = {'qChosenDiff'}; %which value to plot for latencies (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
valType_choice    = {'qDiff'}; %which value to plot for choice (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
latencyType       = {'trialInit_thresh'}; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press)
perfThresh        = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename      = fullfile(whereAreWe('figureCode'), 'processed_data'); % where to save data 
cutoff            = inf; % upper cutoff for latency
lowCutoff         = 0; % lower cutoff for latency 
qFile             = 'qLearn_session_all.mat'; %name of file with q-values
fext              = 'fig1';
weightSpreadsheet = fullfile(whereAreWe('figurecode'), 'raw_data','sessionWeights.xlsx');
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort? 
intervalThresh    = 300; 
%% Which functions to run

extractAllData = 0; % concatenate all data 

extractFlag(1) = 0; % (1)extract latencies and choice in value bins and by previous outcome (Figure 1C,H,B,G); 
plotFlag(1)    = 1; 

extractFlag(2) = 0; %(2)plot Q-model parameters (Figure 1E);
plotFlag(2)    = 1;  


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

%% Extract data

if extractFlag(1)
     ctrl_valueXsex(zscoreFlag,binNum,binNum_choice, cohort_opto,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext,behaviorTable);
end

if extractFlag(2)
    extractQParams(aids_f, aids_m, perfThresh, sessionLength, qFile,ext,basefilename);
end



%% Plot latency x value x sex 

if plotFlag(1)
    
   
    % try to load extracted value, if not saved, extract
    try
        load(fullfile(basefilename,fext,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,fext,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    catch
        ctrl_valueXsex(zscoreFlag,binNum,binNum_choice, cohort_opto,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext);
        load(fullfile(basefilename,fext,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,fext,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    end
    
    stats_latency.posthoc     = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins);
    stats_latency.lme = lme_latency_weightXtrialXsex(aids_f,aids_m,behaviorTable,valType);
    
    try
        load(fullfile(basefilename,fext,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
        load(fullfile(basefilename,fext,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
    catch
        ctrl_valueXsex(zscoreFlag,binNum,binNum_choice, cohort_opto,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext);
        load(fullfile(basefilename,fext,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
        load(fullfile(basefilename,fext,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
    end
    choiceXvalueXsex_plot(choice_f, choice_m, valType_choice,binsChoice);
    stats_choice.lme = lme_choice_weightXtrialXsex(aids_f,aids_m,behaviorTable,valType_choice);
end

if plotFlag(2)
    try 
        load(fullfile(basefilename, fext,sprintf('qParams_perfThres%d_%s_%s_%s.mat', perfThresh,sessionLength,qFile,ext)))
    catch
        extractQParams(aids_f, aids_m, perfThresh, sessionLength, qFile,ext,fullfile(basefilename,fext));
        load(fullfile(basefilename,fext, sprintf('qParams_perfThres%d_%s_%s_%s.mat', perfThresh,sessionLength,qFile,ext)))
    end
    plotQValueParameters(params_f,params_m,qFile);
end

