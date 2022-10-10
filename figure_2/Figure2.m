% Figure 2d-f

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
qFile           = 'qLearn_session_all_2022.mat'; % file with value estimates from Q-learning model
% latency cutoffs
cutoff          = inf; 
lowCutoff       = 0; 
sessionLength   = 'long';
intervalThresh  = 300; 

basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data');




%% Load behavior data

load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');

%% Plot parameters from inverse gaussian distribution for laser and non-laser trials (figure 2d and S7a)

% Parameters
groupingVar     = 'none'; % for fitting inverse gaussian distribution
ver             = 1; % which parameterization of the inverse gaussian

% Figure 2d 
laserType       = 1;
load(fullfile(basefilename,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))

stats_invGauss_outcome=plotInverseGaussian_parameters(fits,{'f';'m'});
plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});


%% Plot mean latencies in value quantile bins
% Figure 2e, SuppFig 7b
laserType      = 1;
for nl = 1:numel(latencyType)
    for nv = 1:numel(valType)
        stats_latencyQuant = latencyQuant_opto(cohort,behaviorTable,valType{nv},binNum,laserType,zscoreFlag,latencyType{nl});
        anovaIn = stats_latencyQuant.glme_diff_anova;
        coeffIn = stats_latencyQuant.glme_diff_coeff;
        tbl = statsTable(anovaIn,coeffIn,0,1,fullfile(basefilename,'stats_tables','ST7.csv'));
    end
end



%% Plot prob(choice = right) in value quantile bins
% Figure 2F, SuppFig 7C
laserType     = 1;
for nv = 1:numel(valType_choice)
    stats_choiceQuant = choice_opto(behaviorTable,valType_choice{nv},binNum_choice,laserType);
    anovaIn = stats_choiceQuant.glme_diff_anova;
    coeffIn = stats_choiceQuant.glme_diff_coeff;
    tbl = statsTable(anovaIn,coeffIn,0,1,fullfile(basefilename,'stats_tables','ST8.csv'));
end

%% Plot laser effect in first and second half of the session (figure S10)
latencyType     = {'trialInit'}; % trialInit (trial initiation latency), leverPress (lever press latency)

for nl = 1:numel(latencyType)
    for nv = 1:numel(valType)
        stats_time=latencyQuant_optoXtime(cohort,behaviorTable,valType{nv},binNum,laserType,zscoreFlag,latencyType{nl});
    end
end

%% Plot control session value modulation v. effect of inhibtion
% Figure S9
% Parameters
groupingVar     = 'none'; % for fitting inverse gaussian distribution
ver             = 1; % which parameterization of the inverse gaussian
[corr_coeff,p_corr]=figureS9a(behaviorTable,cohort,ext,zscoreFlag,groupingVar,laserType);