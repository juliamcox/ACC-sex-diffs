function fit = extractPredictors(fit,tdtEvents,valueFname)
fbasename_bs = fullfile(whereAreWe('data'),'basis_sets'); % location of the basis set functions

% load value predictors and data partition
load(valueFname); 
%% Create possible predictors 
choice = tdtEvents.ipsiChoice;
reward = tdtEvents.reward;
stay = [NaN tdtEvents.choice(1:end-1)==tdtEvents.choice(2:end)]';
stay(stay==0) = -1;
stay(find(tdtEvents.choice==-1)+1) = NaN;
stay(find(tdtEvents.choice==-1)) = NaN;
prevChoice = [NaN choice(1:end-1)];

% qDiff into ipsi-contra (rather than right - left) 
if sum(tdtEvents.ipsiChoice==tdtEvents.choice) == numel(tdtEvents.choice) 
    qDiff = qLeft-qRight;
    qIpsi = qLeft;
    qContra = qRight;
elseif sum(tdtEvents.ipsiChoice==tdtEvents.choice)==0 || sum(tdtEvents.ipsiChoice==tdtEvents.choice)==sum(tdtEvents.choice==-1)
    qDiff = qRight-qLeft;
    qIpsi = qRight;
    qContra = qLeft;
else
    keyboard
end
qDiff_rl = qRight - qLeft;
latency = tdtEvents.nosePokeEntryTimestamps-tdtEvents.trialStartTimestamps; % trial initation latency 
if size(latency,2)>size(latency,1)
    latency=latency';
end

latencyBins = prctile(latency,linspace(0,100,4));
[~,~,latencyQuant] = histcounts(latency,latencyBins);
latency = zscore(latency);

ntrials = numel(tdtEvents.choice);
fit.dataPartition = dataPartition;

% find indices of multiple rewarded and unrewarded trials, if needed

x=1;
multiFlag = 1;
idx = find(reward==1);
while multiFlag == 1
    if ~isempty(idx)
        d = [0 diff(idx)];
        idx2 = idx(d~=1);
        if ~isempty(idx2)
            multiRewIdx{x} = idx2;
            x = x+1;
            idx = idx(~ismember(idx,idx2));
        else
            multiFlag = 0;
        end
        clear d
    else
        multiFlag = 0;
    end
end
x=1;
multiFlag = 1;
idx = find(reward==0);
while multiFlag == 1
    if ~isempty(idx)
        d = [0 diff(idx)];
        idx2 = idx(d~=1);
        if ~isempty(idx2)
            multiNrewIdx{x} = idx2;
            x = x+1;
            idx = idx(~ismember(idx,idx2));
        else
            multiFlag = 0;
        end
        clear d
    else
        multiFlag = 0;
    end
end

multiRew = zeros(size(tdtEvents.choice))';
multiNRew = zeros(size(tdtEvents.choice))';
for nc = 1:numel(multiRewIdx)
    multiRew(multiRewIdx{nc}) = nc;
end
for nc = 1:numel(multiNrewIdx)
    multiNRew(multiNrewIdx{nc}) = nc;
end
    
prevReward = [NaN tdtEvents.reward(1:end-1)];    

%% Create predictor vectors 
for ne = 1:numel(fit.eventNames)
    % Event predictors 
    fit.event_times{ne} = nan(ntrials,1); % create ntrial x 1 vector of NaNs 
    
    try
        thisEvent = eval(sprintf('tdtEvents.%s',fit.eventNames{ne}));
    catch
        if contains(fit.eventNames{ne},'leverPressIpsi')
            thisEvent = tdtEvents.leverPress(tdtEvents.ipsiChoice==1);
        elseif contains(fit.eventNames{ne},'leverPressContra')
            thisEvent = tdtEvents.leverPress(tdtEvents.ipsiChoice==0);
        elseif contains(fit.eventNames{ne},'nosePokeEntryRew')
            thisEvent = tdtEvents.nosePokeEntry(prevReward==1);
        elseif contains(fit.eventNames{ne},'nosePokeEntryNRew')
            thisEvent = tdtEvents.nosePokeEntry(prevReward==0&prevChoice~=-1);
        else
            keyboard
        end
    end
    
    if contains(fit.eventNames{ne},'Ipsi')
        fit.event_times{ne}(tdtEvents.ipsiChoice==1) = thisEvent;
    elseif contains(fit.eventNames{ne},'Contra')
        fit.event_times{ne}(tdtEvents.ipsiChoice==0) = thisEvent;
    elseif contains(fit.eventNames{ne},'CSPlus')
        fit.event_times{ne}(tdtEvents.reward==1) = thisEvent;
    elseif contains(fit.eventNames{ne},'CSMinus')
        fit.event_times{ne}(tdtEvents.reward==0&tdtEvents.choice~=-1) = thisEvent;
    elseif contains(fit.eventNames{ne}, 'nosePokeEntryRew')
        fit.event_times{ne}(prevReward==1) = thisEvent;
    elseif contains(fit.eventNames{ne},'nosePokeEntryNRew')
        fit.event_times{ne}(prevReward==0&prevChoice~=-1) = thisEvent;
    elseif numel(eval(sprintf('tdtEvents.%s',fit.eventNames{ne}))) == ntrials
        fit.event_times{ne} = thisEvent;
    else
        keyboard
    end
    
    fit.event_times{ne}(tdtEvents.choice==-1) = NaN; % remove omitted trials
   
    % Gain predictors 
    for ng = 1:numel(fit.gainNames{ne})
        try fit.gain_preds{ne,ng} = eval(fit.gainNames{ne}{ng});
        catch
            if contains(fit.gainNames{ne}{ng},'next')
                if contains(fit.gainNames{ne}{ng},'_sq')
                    fit.gain_preds{ne,ng} = [eval(sprintf('%s(2:end).^2;',fit.gainNames{ne}{ng}(1:end-9))); NaN];
                else
                    fit.gain_preds{ne,ng} = [eval(sprintf('%s(2:end)',fit.gainNames{ne}{ng}(1:end-5))); NaN];
                end
            elseif contains(fit.gainNames{ne}{ng},'prev')
                if contains(fit.gainNames{ne}{ng},'_sq')
                    fit.gain_preds{ne,ng} = [NaN; eval(sprintf('%s(1:end-1).^2;',fit.gainNames{ne}{ng}(1:end-9)))];
                else
                    fit.gain_preds{ne,ng} = [NaN; eval(sprintf('%s(1:end-1)',fit.gainNames{ne}{ng}(1:end-5)))];
                end
            elseif contains(fit.gainNames{ne}{ng},'_sq')
                fit.gain_preds{ne,ng} = eval(sprintf('%s.^2;',fit.gainNames{ne}{ng}(1:end-4)));
            else
                keyboard
            end
        end   
    end
    fit.gain_preds{ne,ng+1} = ones(size(fit.gain_preds{ne,ng})); % add intercept
end

%% Load/make basis set for each event if necessary
for ne = 1:numel(fit.eventNames)
    
    try load(fullfile(fbasename_bs, sprintf('bs_%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate)))));
        bs{ne} = eval(sprintf('bs_%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate))));
    catch
        range = [-fit.eventDur(ne,1)*fit.frameRate,fit.eventDur(ne,2)*fit.frameRate];
        order = 4; %order of bsplines (cubic)
        basisobj = create_bspline_basis(range, fit.nbasis(ne), order); %create basis function
        eval(sprintf('bs_%s =  getbasismatrix(range(1):range(2), basisobj);', sprintf('%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate))))); %extract matrix of basis set
        save(fullfile(fbasename_bs, sprintf('bs_%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate)))), sprintf('bs_%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate))));
        bs{ne} = eval(sprintf('bs_%dsec_%d_%s',fit.eventDur(ne,1)+fit.eventDur(ne,2),fit.nbasis(ne),num2str(round(fit.frameRate))));
    end
end
fit.bs = bs;

%% Generate trial indices
if sum(contains(fit.eventNames,'nosePokeEntry'))>0
    timeBack = max(fit.eventDur(contains(fit.eventNames,'nosePokeEntry'),1))*fit.frameRate;
else
    timeBack = 2*fit.frameRate;
end
trialIdx = arrayfun(@(x,y) x:y, tdtEvents.nosePokeEntry(1:end-1)-timeBack, tdtEvents.nosePokeEntry(2:end)-1,'UniformOutput',false);
trialIdx = cat(2,trialIdx,{tdtEvents.nosePokeEntry(end)-timeBack:size(dff,1)});
fit.trialIdx = trialIdx;

