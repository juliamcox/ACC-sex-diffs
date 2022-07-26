 function latencyXstaySwitchXsex(zscoreFlag,perfThresh,cohort, cutoff, lowCutoff,aids_m,aids_f,fext,behaviorTable)

%%% Figure 1 %%%


% Extract trial initiation time, lever press latency or nose poke withdrawal as a function of stay v switch and previous outcome 

%%% Inputs
% zscoreFlag: 0 or 1 to zscore resposnse times
% sessionLength: session length: long, short 
% aids_m and aids_f: if don't input animal lists, defaults to combining
% nphr and yfp 
% fext: extension for save location 

%%% Dependencies
% Stan Q model converted to matlab with qValue_individual_fromMat.mat
% generateAnimalList: to create subject lists
% whereAreWe: to create file path
% value_ctrl: embedded function to extract latency information 



%% Parameters 

ipsiFlag = 0;

if nargin < 8  
   % Generate male, female subject  lists
    aids_m   = generateAnimalList([cohort '_male']);
    aids_m   = cat(1,aids_m,generateAnimalList([cohort '_yfp_male']));
    aids_f   = generateAnimalList([cohort '_female']);
    aids_f   = cat(1,aids_f,generateAnimalList([cohort '_yfp_female']));
    fext     = 'fig1';
end
basefilename = fullfile(whereAreWe('figurecode'),'processed_data');
savehere = fullfile(whereAreWe('figurecode'), 'processed_data',fext);
if ~isdir(savehere)
    mkdir(savehere)
end
%% Extract and save latency x value 

[latency_m] = staySwitch_ctrl(aids_m,behaviorTable,zscoreFlag,cutoff,lowCutoff);
save(fullfile(savehere,['ctrlLatencyStaySwitch_m_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_m')

[latency_f] = staySwitch_ctrl(aids_f,behaviorTable,zscoreFlag,cutoff,lowCutoff);
save(fullfile(savehere,['ctrlLatencyStaySwitch_f_zscore' num2str(zscoreFlag) '_' cohort '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_f')

end


%% Extract trial initiation latency X value
function [latency] = staySwitch_ctrl(aids,behaviorTable,zscoreFlag,cutoff,lowCutoff)

%aids: list of animal IDs
%ext: file extension (LAS,UNI et) 
%epochs: laser conditions to extract 
%zscoreFlag: 0 or 1 to zscore resposnse times 
%qFile: which q values to extract
%binNum: binNum for value
%sessionLength: session length


basefilename = fullfile(whereAreWe('bucket'), 'Operant');
latencyTypes = {'trialInit';'leverPress';'trialInit_thresh'}; 

% Initialize variables for response times X value 
for nl = 1:numel(latencyTypes)
    eval(sprintf('latency.%s_stay = zeros(numel(aids),1);',latencyTypes{nl}));
    eval(sprintf('latency.%s_switch = zeros(numel(aids),1);',latencyTypes{nl}));
    
    eval(sprintf('latency.%s_stayRew = zeros(numel(aids),1);',latencyTypes{nl}));
    eval(sprintf('latency.%s_switchRew = zeros(numel(aids),1);',latencyTypes{nl}));
    
    eval(sprintf('latency.%s_stayNRew = zeros(numel(aids),1);',latencyTypes{nl}));
    eval(sprintf('latency.%s_switchNRew = zeros(numel(aids),1);',latencyTypes{nl}));
end



for na = 1:numel(aids)
    
    fprintf('Processing %s...\n', aids{na})
    
    aIdx = find(strcmp(behaviorTable.aID,aids{na}));
    
    
    % remove trials > cutoff
    clear idx
    for nl = 1:numel(latencyTypes)
        eval(sprintf('idx{nl} = find(behaviorTable.%s(aIdx)<=cutoff&behaviorTable.%s(aIdx)>lowCutoff);',latencyTypes{nl},latencyTypes{nl}))
        eval(sprintf('idx_long{nl} = find(behaviorTable.%s(aIdx)>cutoff);',latencyTypes{nl}))
        eval(sprintf('latency.%s_longPer(na) = numel(idx_long{nl})/numel(behaviorTable.%s(aIdx));',latencyTypes{nl},latencyTypes{nl}));
    end
    
    % Divide response times by stay switch, rewarded unrewarded
   % stay = [nan stay(1:end-1)];
    stayIdx = find(behaviorTable.stay(aIdx)==1);
    switchIdx = find(behaviorTable.stay(aIdx)==0);
    stayRewIdx = find(behaviorTable.stay(aIdx)==1&behaviorTable.previousReward(aIdx)==1);
    switchRewIdx = find(behaviorTable.stay(aIdx)==0&behaviorTable.previousReward(aIdx)==1);
    stayNRewIdx = find(behaviorTable.stay(aIdx)==1&behaviorTable.previousReward(aIdx)==0);
    switchNRewIdx = find(behaviorTable.stay(aIdx)==0&behaviorTable.previousReward(aIdx)==0);
    
    
    if zscoreFlag==1
        for nl = 1:numel(latencyTypes)
            thisIdx = stayIdx(ismember(stayIdx,idx{nl}));
            eval(sprintf('latency.%s_stay(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchIdx(ismember(switchIdx,idx{nl}));
            eval(sprintf('latency.%s_switch(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayRewIdx(ismember(stayRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayRew(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchRewIdx(ismember(switchRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchRew(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayNRewIdx(ismember(stayNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayNRew(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchNRewIdx(ismember(switchNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchNRew(na) = nanmean(behaviorTable.%s_zscore(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_trials{na} = behaviorTable.%s_zscore(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
        end
    elseif zscoreFlag == 0
        for nl = 1:numel(latencyTypes)
            thisIdx = stayIdx(ismember(stayIdx,idx{nl}));
            eval(sprintf('latency.%s_stay(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchIdx(ismember(switchIdx,idx{nl}));
            eval(sprintf('latency.%s_switch(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayRewIdx(ismember(stayRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayRew(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchRewIdx(ismember(switchRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchRew(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayNRewIdx(ismember(stayNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayNRew(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchNRewIdx(ismember(switchNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchNRew(na) = nanmean(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_trials{na} = behaviorTable.%s(aIdx(thisIdx));', latencyTypes{nl}, latencyTypes{nl}));
        end
    elseif zscoreFlag == 2
        for nl = 1:numel(latencyTypes)
            thisIdx = stayIdx(ismember(stayIdx,idx{nl}));
            eval(sprintf('latency.%s_stay(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stay_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchIdx(ismember(switchIdx,idx{nl}));
            eval(sprintf('latency.%s_switch(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switch_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayRewIdx(ismember(stayRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayRew(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayRew_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchRewIdx(ismember(switchRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchRew(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchRew_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = stayNRewIdx(ismember(stayNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_stayNRew(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_stayNRew_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
            thisIdx = switchNRewIdx(ismember(switchNRewIdx,idx{nl}));
            eval(sprintf('latency.%s_switchNRew(na) = nanmean(log(behaviorTable.%s(aIdx(thisIdx))));', latencyTypes{nl}, latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_count(na) = numel(thisIdx);', latencyTypes{nl}));
            eval(sprintf('latency.%s_switchNRew_trials{na} = log(behaviorTable.%s(aIdx(thisIdx)));', latencyTypes{nl}, latencyTypes{nl}));
        end
    end
end
    
    
    
  
end