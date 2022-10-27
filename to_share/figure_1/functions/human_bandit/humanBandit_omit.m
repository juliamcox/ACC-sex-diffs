function [initOmit_f,omit_f,initOmit_m,omit_m]=humanBandit_omit(cohort,binNum,zscoreFlag,perfThresh,sideThresh,lowCutoff,highCutoff)

% cohort: which subjects to plot
% qLoc: location of q-learning stan fit 
basefilename = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit');
load(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f');
load(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',binNum,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m');

%% Load each subject's data extract value bins 
sList_f = latency_f.aList; 
sList_m = latency_m.aList;



for ns = 1:numel(sList_f)
    load(fullfile(whereAreWe('figurecode'),'processed_data','human_bandit',[cohort '_f'],sList_f{ns}(8:end))); % load behavior data
    initOmit_f(ns) = mean(data.initOmitIdx);
    omit_f(ns) = mean(data.omitIdx);   
end

for ns = 1:numel(sList_m)
    load(fullfile(whereAreWe('figurecode'),'processed_data','human_bandit',[cohort '_m'],sList_m{ns}(8:end))); % load behavior data
    initOmit_m(ns) = mean(data.initOmitIdx);
    omit_m(ns) = mean(data.omitIdx);   
end
