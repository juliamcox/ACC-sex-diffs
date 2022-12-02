function latencyXtime(behaviorTable,basefilename,binNum,zscoreFlag,perfThresh,qFile,cohort)


behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

%% Load data
load(fullfile(basefilename,['ctrlLatencyXtime_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
load(fullfile(basefilename,['ctrlLatencyXtime_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))

%% Plot
figure();

aids_m = unique(behaviorTable.aID(behaviorTable.female==0));
aids_f = unique(behaviorTable.aID(behaviorTable.female==1));

for nt = 1:numel(latency_m)
    subplot(1,numel(latency_m),nt); hold on
    shadedErrorBar(1:binNum+1,mean(latency_m{nt}.trialInit_qChosenDiff_quant,1),std(latency_m{nt}.trialInit_qChosenDiff_quant)./sqrt(numel(aids_m)),'lineprops',{ 'Color',plotParams.maleC});
    shadedErrorBar(1:binNum+1,mean(latency_f{nt}.trialInit_qChosenDiff_quant,1),std(latency_f{nt}.trialInit_qChosenDiff_quant)./sqrt(numel(aids_f)),'lineprops',{ 'Color',plotParams.femaleC}); 
    if nt == 2
        xlabel('QChosen - QUnchosen _quantile')
    end
    if nt == 1
        ylabel('Trial initiation latency (s)')
    end
end

for nt = 1:numel(latency_m)
    subplot(1,numel(latency_m),nt)
    set(gca,'YLim',[0 13],'XLim',[.5 binNum+.5])
end
end