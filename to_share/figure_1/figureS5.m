
clear all
close all
%% Parameters
cohort          = 'all';
binNum          = 5;
binNum_choice   = 11;
zscoreFlag      = 1; 
valType         = {'QChosenDiff'};
valType_choice  = {'QDiff'};
latencyType     = {'trialStart'}; 
highCutoff      = inf;
lowCutoff       = 0; 
ageLow          = 1;
ageHigh         = 71;
perfThresh      = .1;
sideThresh      = .4;

qFileLoc        = fullfile(whereAreWe('figurecode'),'processed_data','dataForQ','human_bandit','all');
basefilename    = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit');



%% Which functions to run

convertData = 0; % convert and sort data by sex
extractQ    = 0; % compile data for fitting Q-learning model
extractVal  = 0; % extract trial-by-trial value estimates after fitting q-model 
extractLatency = 0; % [extract latencies in value bins, extract choice in value bins]
plotLatency    = 1; % plot latencies (Figure S5f-g)
extractChoice = 0; % extract choice 
plotChoice = 0; % plot choice (Figure S5c-d)
plotQ = 0; % plot parameter estimates from q-model (Figure S5b)
plotLatencyXage = 0; % plot latency vs. age for males and females (Figure S5e)
plotScatter = 1; % plot latency x age scatter plot (Figure S5e)
%% convertData

if convertData
   humanBandit_dataAllocation('cohort1'); 
   humanBandit_dataAllocation('cohort2'); 
end
%% Set up for q-learning model 
if extractQ 
   dataForQ_humanBandit({'cohort1';'cohort2'},qFileLoc); 
end

%% Extract value
if extractVal
   % females
   flist = dir(fullfile(whereAreWe('figurecode'),'processed_data','human_bandit','all_f','*.mat'));
   flist = {flist(:).name};
   flist(contains(flist,'qLearn')) = [];
   extractValue_humanBandit(flist,'all_f',qFileLoc);
   % males
   flist = dir(fullfile(whereAreWe('figurecode'),'processed_data','human_bandit','all_m','*.mat'));
   flist = {flist(:).name};
   flist(contains(flist,'qLearn')) = [];
   extractValue_humanBandit(flist,'all_m',qFileLoc);
end
%% Extract data 

if extractLatency
    latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
end
if extractChoice
    choiceXvalue_humanBandit_sexDiff(cohort,binNum_choice,perfThresh,sideThresh);
end

%% Plot q-values (Figure S5e)
if plotQ 
    try
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'alpha','beta','side','stay');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'alpha_m','beta_m','side_m','stay_m');
    catch
        latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');
    end
    plotQParams(alpha,alpha_m,beta,beta_m,stay,stay_m,side,side_m)
end


%% Plot latency x value x sex 

if plotLatency
    try
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');
    catch
        latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');
    end
   stats = plotLatencyXValue_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag); 
   stats{1}.postHoc_age = humanBandit_valueXage(latency_f,latency_m,valType,latencyType,binNum,zscoreFlag);
  stats_mdl = latencyXValue_linReg_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag);
end

%% Plot choice x value x sex 
if plotChoice
    try
        load(fullfile(basefilename,cohort, sprintf('choiceXvalue_f_bins%d_perfThresh%s_sideThresh%s.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_f','bins');
        load(fullfile(basefilename,cohort, sprintf('choiceXvalue_m_bins%d_perfThresh%s_sideThresh%s.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_m','bins');
    catch
        choiceXvalue_humanBandit_sexDiff(cohort,binNum_choice,perfThresh,sideThresh);
        load(fullfile(basefilename,cohort, sprintf('choiceXvalue_f_bins%d_perfThresh%s_sideThresh%s.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_f','bins');
        load(fullfile(basefilename,cohort, sprintf('choiceXvalue_m_bins%d_perfThresh%s_sideThresh%s.mat',binNum_choice,num2str(perfThresh),num2str(sideThresh))),'choice_m','bins');
    end
    plotChoiceXValue_humanBandit_sexDiff(choice_f,choice_m,bins,valType_choice);
    humanBandit_valueXage_choice(choice_f,choice_m,valType_choice,binNum_choice)
    stats_mdl_choice=choiceXValue_linReg_humanBandit_sexDiff(choice_f,choice_m,bins,valType_choice,zscoreFlag);
end

%% Plot latency x age scatter
zscoreFlag = 0;
if plotScatter
    try
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');
    catch
        latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins');
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins');
    end
    [fit_m,fit_f] = humanBandit_ageXlatency_scatter(latency_f,latency_m);
end