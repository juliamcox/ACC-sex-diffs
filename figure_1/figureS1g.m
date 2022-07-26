    %function figureS1g
%% Parameters
plotFlag        = 1;
extractFlag     = 0;

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

ext             = 'TAB'; % experiment (control sessions)
sessionLength   = 'long';% session length (alternative: 'short';'both')

binNum          = 4; %number of quantile bins for latency
binNum_choice   = 9; %number of quantile bins for choice
zscoreFlag      = 0; %zscore latency?
valType         = {'qChosenDiff'}; %which value to plot (alternative: 'qChosen','qTot','qDiff','qRight','qLeft)
valType_choice  = {'qDiff'};
latencyType     = {'trialInit'}; %which latency to plot (alternative: 'leverPress' (time from lever presentation to lever press), 'withdrawal' (nose poke entry to exit) 
perfThresh      = 0.1; % stay probability reward - stay probability no reward performance threshold 
basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data'); % where to save data 
cutoff          = inf; % upper cutoff for latency
lowCutoff       = 0; % lower cutoff for latency 
qFile           = 'qLearn_session_all_2022.mat'; %name of file with q-values
bandwidth       = .5; % bandwidth for kernel density function
ptiles          = [20:20:100]; % 
fext            = 'fig1';
cohort_opto     = 'ACC_DMS_nphr';
intervalThresh  = 300;
savehere = fullfile(whereAreWe('figurecode'), 'processed_data',fext);


% load all behavior data
load(fullfile(fullfile(whereAreWe('figurecode'), 'processed_data'),sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)));
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);
%%
if extractFlag
    latencyXvalueXsex_time(zscoreFlag,binNum, cohort,sessionLength,perfThresh,qFile,aids_m,aids_f,fext,behaviorTable,flist_f,flist_m)
    load(fullfile(savehere,['ctrlLatencyXtime_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_' qFile]))
    load(fullfile(savehere,['ctrlLatencyXtime_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
else 
    try load(fullfile(savehere,['ctrlLatencyXtime_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
        load(fullfile(savehere,['ctrlLatencyXtime_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
    catch
        latencyXvalueXsex_time(zscoreFlag,binNum, cohort,sessionLength,perfThresh,qFile,aids_m,aids_f,fext,behaviorTable,flist_f,flist_m)
        load(fullfile(savehere,['ctrlLatencyXtime_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
        load(fullfile(savehere,['ctrlLatencyXtime_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]))
    end
end


if plotFlag
    plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
    figure();
    
    for nt = 1:numel(latency_m)
        subplot(1,numel(latency_m),nt); hold on
        shadedErrorBar(1:binNum+1,mean(latency_m{nt}.trialInit_qChosenDiff,1),std(latency_m{nt}.trialInit_qChosenDiff)./sqrt(numel(aids_m)),'lineprops',{ 'Color',plotParams.maleC});
        shadedErrorBar(1:binNum+1,mean(latency_f{nt}.trialInit_qChosenDiff,1),std(latency_f{nt}.trialInit_qChosenDiff)./sqrt(numel(aids_f)),'lineprops',{ 'Color',plotParams.femaleC});
        %
        %     plot(repmat(1:binNum+1,size(latency_f{nt}.trialInit_qChosenDiff_quant,1),1)',latency_f{nt}.trialInit_qChosenDiff_quant','Color',[plotParams.femaleC .25])
        %     plot(repmat(1:binNum+1,size(latency_m{nt}.trialInit_qChosenDiff_quant,1),1)',latency_m{nt}.trialInit_qChosenDiff_quant','Color',[plotParams.maleC .25])
        %
        %     errorbar(1:binNum+1,mean(latency_m{nt}.trialInit_qChosenDiff_quant,1),std(latency_m{nt}.trialInit_qChosenDiff_quant)./sqrt(numel(aids_m)),'Color',plotParams.maleC,'LineWidth',1.5,'CapSize',0)
        %
        %     errorbar(1:binNum+1,mean(latency_f{nt}.trialInit_qChosenDiff_quant,1),std(latency_f{nt}.trialInit_qChosenDiff_quant)./sqrt(numel(aids_f)),'Color',plotParams.femaleC,'LineWidth',1.5,'CapSize',0)
        
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