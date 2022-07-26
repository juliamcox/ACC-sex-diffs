 function ctrl_valueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext,behaviorTable)

%%% Figure 1 %%%


% Extract trial initiation time, lever press latency or nose poke withdrawal as a function of Q-value
% Extract choice as a function of Q-value
% Extract choice and latencies as function of previous outcome


%%% Inputs
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% zscoreFlag: 0 or 1 to zscore resposnse times
% sessionLength: session length: long, short 
% aids_m and aids_f: if don't input animal lists, defaults to combining
% nphr and yfp 
% fext: extension for save location 

%% Parameters 


if nargin < 10  
   % Generate male, female subject  lists
    aids_m   = generateAnimalList([cohort '_male']);
    aids_m   = cat(1,aids_m,generateAnimalList([cohort '_yfp_male']));
    aids_f   = generateAnimalList([cohort '_female']);
    aids_f   = cat(1,aids_f,generateAnimalList([cohort '_yfp_female']));
    fext     = 'fig1';
end


savehere = fullfile(whereAreWe('figurecode'), 'processed_data', fext);
if ~isdir(savehere)
    mkdir(savehere)
end


intervalThresh = 600;


ids = unique(behaviorTable.aID);
behaviorTable.trialInit_thresh = behaviorTable.trialInit; 
for na = 1:numel(ids)
    thisSession = unique(behaviorTable.session(strcmp(behaviorTable.aID,ids{na})));
    for ns = 1:numel(thisSession)
        sessionIdx = find(behaviorTable.session==thisSession(ns)&strcmp(behaviorTable.aID,ids{na})); 
        threshIdx = sessionIdx(find(behaviorTable.trialInit(sessionIdx)>intervalThresh,1,'first'));
        behaviorTable.trialInit_thresh(threshIdx-1:sessionIdx(end)) = NaN;
    end
end

%% Extract and save latency x value 

[latency_m,choice_m,bins,binsChoice] = value_ctrl(aids_m,zscoreFlag,binNum,binNum_choice,behaviorTable,cutoff,lowCutoff);
save(fullfile(savehere,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_m','bins')
save(fullfile(savehere,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']), 'choice_m','binsChoice')

[latency_f,choice_f,bins,binsChoice] = value_ctrl(aids_f,zscoreFlag,binNum,binNum_choice,behaviorTable,cutoff,lowCutoff);
save(fullfile(savehere,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_f','bins')
save(fullfile(savehere,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']), 'choice_f','binsChoice')

end


%% Extract trial initiation latency X value
function [latency,choice,bins,binsChoice] = value_ctrl(aids,zscoreFlag,binNum,binNum_choice,behaviorTable,cutoff,lowCutoff)

%%% Inputs: %%%
%aids: list of animal IDs
%zscoreFlag: 0, no zscore; 1 zscore resposnse times; 2 log of response times 
%qFile: which q values to extract
%binNum: binNum for value
%binNum_choice
%sessionLength: session length ('long';'short';'both')
%perfThresh: reward-no reward stay probability
%cutoff: high cutoff for latencies (in seconds)
%lowCutoff: low cutoff for latencies (in seconds)

%%% Outputs %%%
% latency: data structure with average latencies by value bin for each subject
% choice: data structure with probability choice = right by value bin for each subject
% bins: data structure with bins and quantile bins for each q-value for latency
% binsChoice: data structure with bins and quantile bins for each q-value for choice
% thisFlist: cell array of file names for each subject  

%% Parmeters
% Generate bins for q-values 
bins.qDiff = linspace(-1,1,binNum+1);
bins.qTot = linspace(0,2,binNum+1);
bins.qChosen = linspace(0,1,binNum+1);
bins.qUnchosen = linspace(0,1,binNum+1);
bins.qChosenDiff = linspace(-1,1,binNum+1);

binsChoice.qDiff = linspace(-1,1,binNum_choice+1);
binsChoice.qTot = linspace(0,2,binNum_choice+1);
binsChoice.qChosen = linspace(0,1,binNum_choice+1);
binsChoice.qUnchosen = linspace(0,1,binNum_choice+1);
binsChoice.qChosenDiff = linspace(-1,1,binNum_choice+1);


% Location of the behavior data 
basefilename = fullfile(whereAreWe('behavior')); 
% Value types to extract
valTypes = {'qTot';'qChosenDiff'}; 
valTypes_choice = {'qDiff'}; 
% Latencies to extract 
latencyTypes = {'trialInit';'leverPress';'trialInit_thresh'}; 

%% Initialize variables for response times X value 
for nv = 1:numel(valTypes)
    for nl = 1:numel(latencyTypes)
        eval(sprintf('latency.%s_%s = zeros(numel(aids),numel(bins.%s));',latencyTypes{nl}, valTypes{nv}, valTypes{nv}));
    end
end

for nl = 1:numel(latencyTypes)
    eval(sprintf('latency.%s_longPer = zeros(numel(aids),1);',latencyTypes{nl})); % initialize variable for proportion of long trials
end

%% Extract latency data for each subject 

for na = 1:numel(aids)
    % find indices for mouse na
    aIdx = find(strcmp(behaviorTable.aID,aids{na})&behaviorTable.laserSession==0);
    % generate stay 
    stay = behaviorTable.stay(aIdx);
    % remove session borders
    sessionBreaks = find(diff(behaviorTable.session(aIdx))~=0);
    stay((sessionBreaks+1)) = NaN;
    choice.stay(na) = nanmean(stay);
    
    clear idx
    for nl = 1:numel(latencyTypes)
       eval(sprintf('idx{nl} = find(behaviorTable.%s(aIdx)<=cutoff&behaviorTable.%s(aIdx)>=lowCutoff);',latencyTypes{nl},latencyTypes{nl})) 
    end
    
    fprintf('Processing %s...\n', aids{na})
  
    % bin values
    for nv = 1:numel(valTypes)
       eval(sprintf('[~,~,%s_bins] = histcounts(behaviorTable.%s(aIdx),bins.%s);',valTypes{nv},valTypes{nv},valTypes{nv}));   
    end
    
    for nv = 1:numel(valTypes_choice)
        eval(sprintf('[~,~,%s_bins_choice] = histcounts(behaviorTable.%s(aIdx),bins.%s);',valTypes_choice{nv},valTypes_choice{nv},valTypes_choice{nv}));
    end
    

    % Divide response times by current trial value for laser and non-laser trials
    for nv = 1:numel(valTypes)
        thisbins = eval(sprintf('bins.%s;',valTypes{nv}));
        for nb = 1:numel(thisbins)
            if zscoreFlag==1
                for nl = 1:numel(latencyTypes)
                    thisIdx = eval(sprintf('find(%s_bins==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s(na,nb) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                    eval(sprintf('latency.%s_%s_count(na,nb) = numel(thisIdx);', latencyTypes{nl}, valTypes{nv}));
                    thisIdx = eval(sprintf('find(behaviorTable.%s_quant(aIdx)==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s_quant(na,nb) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                end
            elseif zscoreFlag == 0
                for nl = 1:numel(latencyTypes)
                    thisIdx = eval(sprintf('find(%s_bins==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s(na,nb) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                    eval(sprintf('latency.%s_%s_count(na,nb) = numel(thisIdx);', latencyTypes{nl}, valTypes{nv}));
                    thisIdx = eval(sprintf('find(behaviorTable.%s_quant(aIdx)==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s_quant(na,nb) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                end
            elseif zscoreFlag == 2
                for nl = 1:numel(latencyTypes)
                    thisIdx = eval(sprintf('find(%s_bins==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s(na,nb) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                    eval(sprintf('latency.%s_%s_count(na,nb) = numel(thisIdx);', latencyTypes{nl}, valTypes{nv}));
                    thisIdx = eval(sprintf('find(behaviorTable.%s_quant(aIdx)==nb)',valTypes{nv}));
                    thisIdx = thisIdx(ismember(thisIdx,idx{nl}));
                    eval(sprintf('latency.%s_%s_quant(na,nb) = nanmean(log(behaviorTable.%s_zscore(aIdx(thisIdx))));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                    
                end
            end
        end
    end
    
    % Extract latency based on previous outcome
    for nl = 1:numel(latencyTypes)
        if zscoreFlag == 1
            thisIdx_rew = find(behaviorTable.previousReward(aIdx)==1);
            thisIdx_nrew = find(behaviorTable.previousReward(aIdx)==0);
            thisIdx_rew = thisIdx_rew(ismember(thisIdx_rew,idx{nl}));
            thisIdx_nrew = thisIdx_nrew(ismember(thisIdx_nrew,idx{nl}));
            eval(sprintf('latency.%s_prevReward(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx_rew)));',latencyTypes{nl},latencyTypes{nl}));
            eval(sprintf('latency.%s_prevNoReward(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx_nrew)));',latencyTypes{nl},latencyTypes{nl}));
        elseif zscoreFlag == 0
            thisIdx_rew = find(behaviorTable.previousReward(aIdx)==1);
            thisIdx_nrew = find(behaviorTable.previousReward(aIdx)==0);
            thisIdx_rew = thisIdx_rew(ismember(thisIdx_rew,idx{nl}));
            thisIdx_nrew = thisIdx_nrew(ismember(thisIdx_nrew,idx{nl}));
            eval(sprintf('latency.%s_prevReward(na) = nanmean(behaviorTable.%s(aIdx(thisIdx_rew)));',latencyTypes{nl},latencyTypes{nl}));
            eval(sprintf('latency.%s_prevNoReward(na) = nanmean(behaviorTable.%s(aIdx(thisIdx_nrew)));',latencyTypes{nl},latencyTypes{nl}));
        elseif zscoreFlag == 2
            thisIdx_rew = find(behaviorTable.previousReward(aIdx)==1);
            thisIdx_nrew = find(behaviorTable.previousReward(aIdx)==0);
            thisIdx_rew = thisIdx_rew(ismember(thisIdx_rew,idx{nl}));
            thisIdx_nrew = thisIdx_nrew(ismember(thisIdx_nrew,idx{nl}));
            eval(sprintf('latency.%s_prevReward(na) = nanmean(log(behaviorTable.%s_zscore(aIdx(thisIdx_rew))));',latencyTypes{nl},latencyTypes{nl}));
            eval(sprintf('latency.%s_prevNoReward(na) = nanmean(log(behaviorTable.%s_zscore(aIdx(thisIdx_nrew))));',latencyTypes{nl},latencyTypes{nl}));
        end
    end
    
    
    %% Divide choice by previous outcome and q-value
    for nv = 1:numel(valTypes_choice)
        thisbins = eval(sprintf('binsChoice.%s;',valTypes_choice{nv}));
        for nb = 1:numel(thisbins)
            % non-laser trials
            thisIdx = eval(sprintf('find(behaviorTable.%s_quant_choice(aIdx)==nb&~isnan(behaviorTable.trialInit_thresh(aIdx)))',valTypes_choice{nv}));
            eval(sprintf('choice.stayProb_%s_quant(na,nb) = nanmean(stay(thisIdx));', valTypes_choice{nv}));
            eval(sprintf('choice.choice_%s_quant(na,nb) = nansum(behaviorTable.choice(aIdx(thisIdx))==0)/numel(thisIdx);', valTypes_choice{nv}));

            eval(sprintf('choice.stayProb_%s(na,nb) = nanmean(stay(%s_bins_choice==nb));', valTypes_choice{nv},valTypes_choice{nv}));
            eval(sprintf('choice.choice_%s(na,nb) = nansum(behaviorTable.choice(aIdx(%s_bins_choice==nb))==0)/sum(%s_bins_choice==nb);', valTypes_choice{nv},valTypes_choice{nv},valTypes_choice{nv}));
        end
        choice.stayProb_prevReward(na) = nanmean(stay(behaviorTable.previousReward(aIdx)==1));
        choice.stayProb_prevNoReward(na) = nanmean(stay(behaviorTable.previousReward(aIdx)==0));
    end
    
    
    
end
end