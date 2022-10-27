
clear all
close all
%% Parameters
cohort          = 'cohort2';
qLoc            = 'dataForQ_cohort2/05-Feb-2022';
% cohort            = 'cohort1';
% qLoc              = 'dataForQ_cohort1/19-Jan-2021';
cohort          = 'all';
qLoc            = 'dataForQ_all/06-Feb-2022';

% 
binNum         =5;
binNum_choice  = 11;
zscoreFlag      =1; 
valType         = {'QChosenDiff'};
valType_choice  = {'QDiff'};
latencyType     = {'trialStart'}; 
highCutoff      = inf;%10000; 
lowCutoff       = 0; 
perfThresh      = 0.1;
sideThresh      = .4;
ageLow          = 1;
ageHigh         = 71;
basefilename    = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit'); 
%% Which functions to run

extractFlag = [1 0]; % [extract latencies in value bins (latency and choice)]
plotFlag    = [1 0]; % plot 

%% Extract value

if extractFlag(1)
    latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
end
if extractFlag(2)
    choiceXvalue_humanBandit_sexDiff(cohort,binNum_choice,ageLow,ageHigh,perfThresh,sideThresh)
end


%% Plot latency x value x sex 

if plotFlag(1)
    try
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)));
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)));
    catch
        latencyXvalue_humanBandit_sexDiff(cohort,binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)));
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)));
    end
   stats= plotLatencyXValue_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag); 
   stats_mdl = latencyXValue_linReg_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag); 
end


if plotFlag(2)
    try
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d.mat',binNum_choice)));
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d.mat',binNum_choice)));
    catch
        choiceXvalue_humanBandit_sexDiff(cohort,qLoc,binNum_choice)
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d.mat',binNum_choice)));
        load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d.mat',binNum_choice)));
    end
    plotChoiceXValue_humanBandit_sexDiff(choice_f,choice_m,bins,valType_choice);
    stats_mdl=choiceXValue_linReg_humanBandit_sexDiff(choice_f,choice_m,bins,valType_choice,zscoreFlag);
end


%% Plot mean latency x age

meanLatency_f = cellfun(@(x) nanmean(x),latency_f.trialStart_trials);
meanLatency_m = cellfun(@(x) nanmean(x),latency_m.trialStart_trials);
plotParams = load(fullfile(whereAreWe('bucket'),'Manuscript_figures','plotParams.mat'));
figure(); hold on
scatter(latency_f.age, meanLatency_f,20,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5)
scatter(latency_m.age, meanLatency_m,20,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5)
l = lsline;
l(1).Color = plotParams.maleC;
l(2).Color = plotParams.femaleC;
l(1).LineWidth = 1;
l(2).LineWidth = 1;
set(gca,'XLim',[18 71])

plotQParams(alpha,alpha_m,beta,beta_m,stay,stay_m,side,side_m)
