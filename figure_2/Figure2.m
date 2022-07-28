% Figure 2 D-F, Supp Fig 5, Supp Fig 6 

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

fext            = 'fig2';


%% Load behavior data
try
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
catch
    weightSpreadsheet = fullfile(whereAreWe('figurecode'), 'raw_data','sessionWeights.xlsx');
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
    extractData(aids_m,aids_f, perfThresh,qFile,weightSpreadsheet,ext,binNum,binNum_choice,cohort,sessionLength,intervalThresh);
    load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');
end

%% Fit inverse gaussian distribution for laser and non-laser trials 

% Parameters
groupingVar     = 'none'; % for fitting inverse gaussian distribution
ver             = 1; % which parameterization of the inverse gaussian


% Figure 2D 
laserType       = 1;
try
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
catch
    fitLatencyDistributions_inverseGaussian(cohort,ext,zscoreFlag,fullfile(basefilename,fext),laserType,valType{1},groupingVar,ver,binNum,behaviorTable,'trialInit_thresh');
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
end
stats_invGauss_outcome=plotInverseGaussian_parameters(fits,{'f';'m'});
stats_invGauss_outcome_yfp=plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});


% Supp Figure 8
laserType       = 2;
try
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
catch
    fitLatencyDistributions_inverseGaussian(cohort,ext,zscoreFlag,fullfile(basefilename,fext),laserType,valType{1},groupingVar,ver,binNum,behaviorTable,'trialInit');
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
end
stats_invGauss_nosePoke=plotInverseGaussian_parameters(fits,{'f';'m'});
plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});


laserType       = 3;
try
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
catch
    fitLatencyDistributions_inverseGaussian(cohort,ext,zscoreFlag,fullfile(basefilename,fext),laserType,valType{1},groupingVar,ver,binNum,behaviorTable,'trialInit');
    load(fullfile(basefilename,fext,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
end
stats_invGauss_ITI=plotInverseGaussian_parameters(fits,{'f';'m'});
plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});


%% Plot mean latencies in value quantile bins
% Figure 1E, SuppFig 5B
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
% Figure 1F, SuppFig 5C
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