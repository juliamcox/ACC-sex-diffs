
close all
clear all
%% Parameters
aids_f = [];
aids_m = [];
aids_f = generateAnimalList('ACC_DMS_nphr_female');
aids_f = cat(1,aids_f,generateAnimalList('ACC_DMS_nphr_yfp_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d1_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d2_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_yfp_female')); 
%aids_f = cat(1,aids_f,generateAnimalList('imaging_female'));

aids_m = generateAnimalList('ACC_DMS_nphr_male');
aids_m = cat(1,aids_m,generateAnimalList('ACC_DMS_nphr_yfp_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d1_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d2_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_yfp_male')); 
%aids_m = cat(1,aids_m,generateAnimalList('imaging_male'));

cohort = 'all';

ext             = 'TAB'; 
sessionLength   = 'long';

binNum          = 4;
binNum_choice   = 9;
zscoreFlag      = 0; 
valType         = {'qChosenDiff'};
valType_choice  = {'qDiff'};
latencyType     = {'trialInit_thresh'}; 
perfThresh      = 0.1; 
basefilename    = fullfile(whereAreWe('figurecode'), 'processed_data');
savehere        = fullfile(whereAreWe('figurecode'), 'processed_data','fig1');

cutoff          = inf; 
lowCutoff       = 0; 
qFile           = 'qLearn_session_all_2022.mat';
cohort_opto     = 'ACC_DMS_nphr';
fext            = 'fig1';
intervalThresh = 300; 

extractFlag = 1; 

try
    % load all behavior data
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
catch
    extractData(aids_m,aids_f, perfThresh,qFile,weightSpreadsheet,ext_laser,binNum,binNum_choice,cohort_opto,sessionLength,intervalThresh);
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
end


behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);


%% extract latencies by previous outcome and stay v. switch
if extractFlag
    latencyXstaySwitchXsex(zscoreFlag,perfThresh,cohort, cutoff, lowCutoff,aids_m,aids_f,fext,behaviorTable);
end
%% load latencies by previous outcome and stay v switch
try
    load(fullfile(savehere,['ctrlLatencyStaySwitch_f_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));
    load(fullfile(savehere,['ctrlLatencyStaySwitch_m_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));
catch
    latencyXstaySwitchXsex(zscoreFlag,sessionLength,perfThresh,cohort, cutoff, lowCutoff,qFile,aids_m,aids_f,fext,cohort_opto,binNum,binNum_choice);
    load(fullfile(savehere,['ctrlLatencyStaySwitch_f_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));
    load(fullfile(savehere,['ctrlLatencyStaySwitch_m_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']));
end

%% Plot stay v switch by previous outcome 
stats=plotStaySwitch_ctrl(latency_f,latency_m,latencyType);

