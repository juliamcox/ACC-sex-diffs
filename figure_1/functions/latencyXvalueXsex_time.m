 function latencyXvalueXsex_time(zscoreFlag,binNum, cohort,sessionLength,perfThresh,qFile,aids_m,aids_f,fext,behaviorTable,flist_f,flist_m)

%%% Figure 1 %%%


% Extract trial initiation time, lever press latency or nose poke withdrawal as a function of Q-value and time in session
% Extract choice and latency as a function of Q-value

%%% Inputs
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% zscoreFlag: 0 or 1 to zscore resposnse times
% sessionLength: session length: long, short 
% aids_m and aids_f: if don't input animal lists, defaults to combining
% nphr and yfp 
% fext: extension for save location 


%% Parameters 


savehere = fullfile(whereAreWe('figurecode'), 'processed_data');
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);


%% Extract and save latency x value 

[latency_m,bins,flist_m,rewardBins_m] = value_ctrl(aids_m,zscoreFlag,binNum,sessionLength,behaviorTable,flist_m);
keyboard
save(fullfile(savehere,fext,['ctrlLatencyXtime_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_' qFile]), 'latency_m','bins','flist_m','rewardBins_m')

[latency_f,bins,flist_f,rewardBins_f] = value_ctrl(aids_f,zscoreFlag,binNum,sessionLength,behaviorTable,flist_f);
save(fullfile(savehere,fext,['ctrlLatencyXtime_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh)  '_' qFile]), 'latency_f','bins','flist_f','rewardBins_f')


end


%% Extract trial initiation latency X value for each time bin
function [latency,bins,thisFlist,rewardBins] = value_ctrl(aids,zscoreFlag,binNum,sessionLength,behaviorTable,flist)
%% inputs
%aids: list of animal IDs
%ext: file extension (LAS,UNI et)
%epochs: laser conditions to extract
%zscoreFlag: 0 or 1 to zscore resposnse times
%qFile: which q values to extract
%binNum: binNum for value
%sessionLength: session length

%% Parameter
basefilename = fullfile(whereAreWe('behavior'));
valTypes = {'qDiff';'qTot';'qChosenDiff'};% Value types to extract 
latencyTypes = {'trialInit';'leverPress'};% Latency types to extract
ntimeBins = 4; %how many time bins
ntimeBins_fine = 30; 

% value bins 
bins.qDiff = linspace(-1,1,binNum+1);
bins.qTot = linspace(0,2,binNum+1);
bins.qChosenDiff = linspace(-1,1,binNum+1);

%% time bins
if contains(sessionLength, 'long')
    rangeStart = 0:floor(7200/ntimeBins):7200;
    if rangeStart(end) ~= 7200
        rangeStart(end+1) = 7200;
    end
    rangeStart_fine = 0:floor(7200/ntimeBins_fine):7200;
    if rangeStart_fine(end) ~= 7200
        rangeStart_fine(end+1) = 7200;
    end
elseif contains(sessionLength,'short')
    rangeStart = 0:floor(3600/ntimeBins):3600;
    if rangeStart(end) ~= 3600
        rangeStart(end+1) = 3600;
    end
    rangeStart_fine = 0:floor(3600/ntimeBins_fine):3600;
    if rangeStart_fine(end) ~= 3600
        rangeStart_fine(end+1) = 3600;
    end
else
    keyboard
end

%% Extract mean latency in value bins for each block of time in session 

% Initialize variables for response times X value
for nt = 1:ntimeBins
    for nv = 1:numel(valTypes)
        for nl = 1:numel(latencyTypes)
            eval(sprintf('latency{nt}.%s_%s = zeros(numel(aids),numel(bins.%s));',latencyTypes{nl}, valTypes{nv}, valTypes{nv}));
            eval(sprintf('latency{nt}.%s_%s_quant = zeros(numel(aids),numel(bins.%s));',latencyTypes{nl}, valTypes{nv}, valTypes{nv}));
            eval(sprintf('latency{nt}.%s_%s_count = zeros(numel(aids),numel(bins.%s));',latencyTypes{nl}, valTypes{nv}, valTypes{nv}));
        end
    end
end

for na = 1:numel(aids)   
    fprintf('Processing %s...\n', aids{na})
    aIdx = find(strcmp(behaviorTable.aID,aids{na}));

    
    % Bin values
    [~,~,qDiff_bins] = histcounts(behaviorTable.qDiff(aIdx),bins.qDiff);
    [~,~,qTot_bins] = histcounts(behaviorTable.qTot(aIdx),bins.qTot);
    [~,~,qChosenDiff_bins] = histcounts(behaviorTable.qChosenDiff(aIdx),bins.qChosenDiff);
    
    
    % extract data for each sub-session
    thisLatency = cell(ntimeBins,1);
    for nt = 1:ntimeBins
        for nl = 1:numel(latencyTypes)
            eval(sprintf('thisLatency{nt}.%s = [];', latencyTypes{nl}));
        end
        for nv = 1:numel(valTypes)
            eval(sprintf('thisLatency{nt}.%s_quant = [];', valTypes{nv}));
        end
        thisLatency{nt}.prevOutcome = [];
        thisLatency{nt}.reward = [];
    end
    session = behaviorTable.session(aIdx);
    sessionIdx = unique(session,'rows','stable');
    
    thisFlist = flist{na}(sessionIdx);
    for nf = 1:numel(thisFlist)
        thisTrials = find(session == sessionIdx(nf));
        thisTable = behaviorTable(strcmp(behaviorTable.aID,aids{na})&behaviorTable.session == sessionIdx(nf),:);
        load(fullfile(whereAreWe('bucket'),'Operant',aids{na}, sprintf('%s_TAB.mat',thisFlist{nf}(2:end-4))));
        data.rt.trialInit = data.rt.nosePoke;
        % indices for value variables (which do not include omitted trials, or first trial) during time bin nt
        %trialIndex_val = arrayfun(@(x,y) find(data.trialStart(data.trialIdx_all)>=x&data.trialStart(data.trialIdx_all)<y), rangeStart(1:end-1), rangeStart(2:end),'UniformOutput',false);
        % indices for data structure variables during time bin nt
        trialIndex = arrayfun(@(x,y) find(data.trialStart>=x&data.trialStart<y), rangeStart(1:end-1), rangeStart(2:end),'UniformOutput',false);
        trialIndex = cellfun(@(x) x(ismember(x,data.trialIdx_all)), trialIndex, 'UniformOutput', false);
        
        trialIndex_fine = arrayfun(@(x,y) find(data.trialStart>=x&data.trialStart<y), rangeStart_fine(1:end-1), rangeStart_fine(2:end),'UniformOutput',false);
        trialIndex_fine = cellfun(@(x) x(ismember(x,data.trialIdx_all)), trialIndex_fine, 'UniformOutput', false);
        
        thisRewardBins(nf,:) = cellfun(@(x) sum(data.reward(x)), trialIndex_fine);
        
        for nt = 1:(ntimeBins)
            for nl = 1:numel(latencyTypes)
                %eval(sprintf('thisLatency{nt}.%s = cat(1,thisLatency{nt}.%s, data.rt.%s(trialIndex{nt})'');', latencyTypes{nl},latencyTypes{nl},latencyTypes{nl}));
                eval(sprintf('thisLatency{nt}.%s = cat(1,thisLatency{nt}.%s, thisTable.%s(ismember(thisTable.trial,trialIndex{nt})));',latencyTypes{nl},latencyTypes{nl},latencyTypes{nl}));
            end
            try
            for nv = 1:numel(valTypes)
                eval(sprintf('thisLatency{nt}.%s_quant = cat(1,thisLatency{nt}.%s_quant, thisTable.%s_quant(ismember(thisTable.trial,trialIndex{nt})));',valTypes{nv},valTypes{nv},valTypes{nv}));
%                 %s_bins(thisTrials(trialIndex{nt})));', valTypes{nv},valTypes{nv},valTypes{nv}));
%                 eval(sprintf('thisLatency{nt}.%s_quant = cat(1,thisLatency{nt}.%s_quant, behaviorTable.%s_quant(thisTrials(trialIndex{nt})));', valTypes{nv},valTypes{nv},valTypes{nv}));
            end
            catch
                keyboard
            end
            thisLatency{nt}.prevOutcome = cat(1,thisLatency{nt}.prevOutcome, data.reward(trialIndex{nt}-1)');
            thisLatency{nt}.reward = cat(1, thisLatency{nt}.reward, data.reward(trialIndex{nt})');
        end
    end
    
    % Divide trials by value
    for nt = 1:ntimeBins
        for nv = 1:numel(valTypes)
            thisbins = eval(sprintf('bins.%s;',valTypes{nv}));
            for nb = 1:numel(thisbins)
                if zscoreFlag==1
                    keyboard
                elseif zscoreFlag == 0
                    for nl = 1:numel(latencyTypes)
                        thisIdx = eval(sprintf('find(thisLatency{nt}.%s_quant==nb)',valTypes{nv}));
                        eval(sprintf('latency{nt}.%s_%s_quant(na,nb) = nanmean(thisLatency{nt}.%s(thisIdx));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                        eval(sprintf('latency{nt}.%s_%s_count(na,nb) = numel(thisIdx);', latencyTypes{nl}, valTypes{nv}));
                    end
                elseif zscoreFlag == 2
                    for nl = 1:numel(latencyTypes)
                        thisIdx = eval(sprintf('find(thisLatency{nt}.%s==nb)',valTypes{nv})); % find trials
                        eval(sprintf('latency{nt}.%s_%s(na,nb) = nanmean(log(thisLatency{nt}.%s(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                        eval(sprintf('latency{nt}.%s_%s_count(na,nb) = numel(thisIdx);', latencyTypes{nl}, valTypes{nv}));
                        thisIdx = eval(sprintf('find(thisLatency{nt}.%s_quant==nb)',valTypes{nv}));
                        eval(sprintf('latency{nt}.%s_%s_quant(na,nb) = nanmean(log(thisLatency{nt}.%s(thisIdx)));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                    end
                end
            end
        end
        latency{nt}.reward_rate(na) = nanmean(thisLatency{nt}.reward);
        latency{nt}.reward_count(na) = nansum(thisLatency{nt}.reward);
        latency{nt}.trial_count(na) = numel(thisLatency{nt}.reward); 
        latency{nt}.prevOutcome(na) = nanmean(thisLatency{nt}.prevOutcome); 
        
      
        
    end
    
    rewardBins{na} = thisRewardBins;
    clear thisRewardBins
    
    
end

end

