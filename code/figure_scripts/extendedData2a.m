

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

latencyType       = {'trialInit_thresh'}; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press)
perfThresh        = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename      = fullfile(whereAreWe('data')); % where to save data 
cutoff            = inf; % upper cutoff for latency
lowCutoff         = 0; % lower cutoff for latency 
qFile             = 'qLearn_session_all.mat'; %name of file with q-values
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort? 
intervalThresh    = 300; % cut session when ITI longer than intervalThresh

extractFlag       = 0; 

%% Load data 


load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);

load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))



%% Plot latency x value x sex 

runReg =0; % flag to fit mixed-effects model 
stats_latency = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins,behaviorTable,aids_f,aids_m,runReg);
if runReg
    anovaIn             = stats_latency.mdl{1}.mdl_anova;
    coefIn              = stats_latency.mdl{1}.mdl_coeff;
    tbl_latency_value   = statsTable(anovaIn, coefIn, 1,0,fullfile(basefilename,'stats_tables','FigS3.csv'));
end

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
tbl = statsTable(anovaIn,coeffIn,1,0,fullfile(whereAreWe('data'),'stats_tables','ST6.csv'));
tbl = statsTable(anovaIn,coeffIn,0,1,fullfile(whereAreWe('data'),'stats_tables','ST6_2.csv'));