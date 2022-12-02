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
latencyType       = {'trialInit_thresh'}; %which latency to plot (alternative: 'leverPress'  (time from lever presentation to lever press) or trialInit (latencies without disengagement cutoff)
perfThresh        = 0.1; % stay probability reward - stay probability no reward performance threshold
basefilename      = whereAreWe('data'); % where is the data
cutoff            = inf; % upper cutoff for latency
lowCutoff         = 0; % lower cutoff for latency
qFile             = 'qLearn_session_all.mat'; %name of file with q-values
fext              = 'fig1';
ext_laser         = 'LAS';
cohort_opto       = 'ACC_DMS_nphr'; % which opto cohort?
intervalThresh    = 300;


%% Load all behavior data and data summaries

load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);
load(fullfile(basefilename,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
load(fullfile(basefilename,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort_opto '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))

%% Plot figure 1b
stats_choice.outcome = choiceXoutcomeXsex_plot(choice_f, choice_m, behaviorTable,aids_f,aids_m);
tbl_choice_outcome   = statsTable(stats_choice.outcome.mdl_anova, stats_choice.outcome.mdl_coeff, 1, 0,fullfile(basefilename,'stats_tables','ST1.csv'));
statsTable(stats_choice.outcome.mdl_anova, stats_choice.outcome.mdl_coeff, 0, 1,fullfile(basefilename,'stats_tables','ST1_tstat.csv'));

%% Plot figure 1c
stats_latency.outcome   = latencyXoutcomeXsex_plot(latency_f, latency_m,latencyType,zscoreFlag,behaviorTable,aids_f,aids_m);

anovaIn = eval(sprintf('stats_latency.outcome.%s.mdl_anova;', latencyType{1}));
coefIn = eval(sprintf('stats_latency.outcome.%s.coeff;', latencyType{1}));
tbl_latency_outcome     = statsTable(anovaIn, coefIn, 1, 0,fullfile(basefilename,'stats_tables','ST2.csv'));
tbl_latency_outcome     = statsTable(anovaIn, coefIn, 0, 2,fullfile(basefilename,'stats_tables','ST2_tstat.csv'));

%% Plot figure 1e
load(fullfile(basefilename ,sprintf('qParams_perfThresh%s_%s_%s_%s', num2str(perfThresh),sessionLength,ext,qFile)))
plotQValueParameters(params_f,params_m);

%% Plot figure 1g

stats_choice = choiceXvalueXsex_plot(choice_f, choice_m, valType_choice,binsChoice,behaviorTable,aids_f,aids_m);
anovaIn             = stats_choice.mdl_anova_qDiff_quant;
coefIn              = stats_choice.mdl_coeff_qDiff_quant;
tbl_choice_value   = statsTable(anovaIn, coefIn, 0, 1,fullfile(basefilename,'stats_tables','ST4.csv'));

%% Plot figure 1 h
runReg = 0; % flag to fit mixed-effects model 
stats_latency.value = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins,behaviorTable,aids_f,aids_m,runReg);
if runReg
    anovaIn             = stats_latency.value.mdl{1}.mdl_anova;
    coefIn              = stats_latency.value.mdl{1}.mdl_coeff;
    tbl_latency_value   = statsTable(anovaIn, coefIn, 0, 1,fullfile(basefilename,'stats_tables','ST5.csv'));
end




