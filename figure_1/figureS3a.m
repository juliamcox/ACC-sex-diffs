

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
valType           = {'qTot'}; %which value to plot for latencies (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
valType_choice    = {'qDiff'}; %which value to plot for latencies (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)

latencyType       = {'trialInit_thresh'}; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press)
perfThresh        = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename      = fullfile(whereAreWe('figureCode'), 'processed_data'); % where to save data 
cutoff            = inf; % upper cutoff for latency
lowCutoff         = 0; % lower cutoff for latency 
qFile             = 'qLearn_session_all_2022.mat'; %name of file with q-values
fext              = 'fig1';
weightSpreadsheet = fullfile(whereAreWe('figurecode'), 'raw_data','sessionWeights.xlsx');
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort? 
intervalThresh    = 300; 

extractFlag       = 0; 

%% Concatenate all data 


try
    % load all behavior data
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
catch
    extractData(aids_m,aids_f, perfThresh,qFile,weightSpreadsheet,ext_laser,binNum,binNum_choice,cohort_opto,sessionLength,intervalThresh);
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
end


behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);

%% Extract data

if extractFlag(1)
     ctrl_valueXsex(zscoreFlag,binNum,binNum_choice, cohort_opto,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext,behaviorTable);
end



%% Plot latency x value x sex 

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

%% Fit model with total value and relative chosen value 

f = sprintf('%s ~ female*qTot + female*qChosenDiff + trial + (1+qTot + qChosenDiff|aID)',latencyType{1}); 
X = behaviorTable;
X.trial = nanzscore(X.trial); 
X.qTot = nanzscore(X.qTot);
X.female = categorical(X.female); 
X.qChosenDiff = nanzscore(X.qChosenDiff);
mdl = fitlme(X,f,'DummyVarCoding','effects');
coeffIn = dataset2cell(mdl.Coefficients);
anovaIn= dataset2cell(anova(mdl,'DFMethod','Satterthwaite'));
tbl = statsTable(anovaIn,coeffIn,1,0,fullfile(whereAreWe('figurecode'),'processed_data','stats_tables','ST6.csv'));