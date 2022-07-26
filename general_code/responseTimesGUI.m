function data = responseTimesGUI(data)
N = length(data.choice);

%%% no laser
%%%% trial type indices

rewIdx = find(data.reward == 1); %reward
nrewIdx = find(data.reward == 0);  %no reward
rewIdx = rewIdx(ismember(rewIdx, data.trialIdx)); %only non-laser trials
nrewIdx = nrewIdx(ismember(nrewIdx, data.trialIdx)); %only non-laser trials

% if rewarded did they choose high probability lever or if unrewarded did they choose the low prob lever?
congIdx = [find(data.reward==1 & data.choice == data.blocks.blockIDAll-1)+1 find(data.reward == 0 & data.choice ~= data.blocks.blockIDAll-1)+1];
% if unrewarded, did they choose the high probability lever or if rewarded did they choose the low prob lever?
incongIdx = [find(data.reward==0 & data.choice == data.blocks.blockIDAll-1)+1 find(data.reward == 1 & data.choice ~= data.blocks.blockIDAll-1)+1];
% non laser trials 
congIdx = congIdx(ismember(congIdx,data.trialIdx));
incongIdx = incongIdx(ismember(incongIdx, data.trialIdx)); 

% Lever press response times by current trial type
rt.leverPress = data.leverPress-data.leverPresentation;
rt.leverPress(data.errors~=0) = NaN;%remove omitted trials 
rt.leverPress(data.leverPress==-10) = NaN;
rt.leverRew = rt.leverPress(rewIdx);
rt.leverNRew = rt.leverPress(nrewIdx);
rt.leverPress_noLaser = rt.leverPress(data.trialIdx);

%%%% Nose poke rt 
rt.nosePoke = data.nosePokeEntry(1:N) - data.trialStart(1:N); 
rt.nosePoke(data.nosePokeEntry==-10) = NaN;
rt.nosePoke_noLaser = data.nosePokeEntry(data.trialIdx)-data.trialStart(data.trialIdx);
%response times by trial type 
rt.noseRew = rt.nosePoke(rewIdx);
rt.noseNRew = rt.nosePoke(nrewIdx);
%%%% Nose poke withdrawal time
rt.withdraw = data.nosePokeExit-data.nosePokeEntry;
rt.withdraw(data.nosePokeExit==-10) = NaN;
rt.withdrawRew = rt.withdraw(rewIdx);
rt.withdrawNRew = rt.withdraw(nrewIdx);
rt.withdraw_noLaser = rt.withdraw(data.trialIdx);

% Response time by previous trial type non-laser trials only
x = (ismember(rewIdx+1, [data.trialIdx_all]));
rt.leverPrevRew = rt.leverPress(rewIdx(x)+1); 
rt.nosePrevRew = rt.nosePoke(rewIdx(x)+1); 
rt.withdrawPrevRew = rt.withdraw(rewIdx(x)+1);
clear x 
x = (ismember(nrewIdx+1, [data.trialIdx_all]));
rt.leverPrevNRew = rt.leverPress(nrewIdx(x)+1); 
rt.nosePrevNRew = rt.nosePoke(nrewIdx(x)+1); 
rt.withdrawPrevNRew = rt.withdraw(nrewIdx(x)+1);

rt.leverPrevCong = rt.leverPress(congIdx);
rt.leverPrevIncong = rt.leverPress(incongIdx);
rt.nosePrevCong = rt.nosePoke(congIdx);
rt.nosePrevIncong = rt.nosePoke(incongIdx);
rt.withdrawPrevCong = rt.withdraw(congIdx);
rt.withdrawPrevIncong = rt.withdraw(incongIdx);


%%% Trial duration
rewIdx = rewIdx(ismember(rewIdx+1, data.trialIdx_all));
nrewIdx = nrewIdx(ismember(nrewIdx+1, data.trialIdx_all));
rt.trialDurRew = data.trialStart(rewIdx+1) - data.trialStart(rewIdx);
rt.trialDurNRew = data.trialStart(nrewIdx+1) - data.trialStart(nrewIdx);
rewIdx = rewIdx(ismember(rewIdx+2, data.trialIdx_all));
nrewIdx = nrewIdx(ismember(nrewIdx+2, data.trialIdx_all));
rt.trialDurPrevRew = data.trialStart(rewIdx+2) - data.trialStart(rewIdx+1);
rt.trialDurPrevNRew = data.trialStart(nrewIdx+2) - data.trialStart(nrewIdx+1);
rt.trialDur = data.trialStart(data.trialIdx(ismember(data.trialIdx+1, data.trialIdx_all))+1) - data.trialStart(data.trialIdx(ismember(data.trialIdx+1, data.trialIdx_all)));

%%% Reward consumption

rt.rewardConsume = data.rewardExit(data.reward==1)-data.rewardEntry(data.reward==1);
rt.rewardConsume_noLaser = data.rewardExit(data.trialIdx(ismember(data.trialIdx, find(data.reward==1))))-data.rewardEntry(data.trialIdx(ismember(data.trialIdx, find(data.reward==1))));

% conditions = {};
% for ni = 1:5
%     conditions = [conditions; ['laserPrevOutcome', num2str(ni)]];
%     conditions = [conditions;['laserOutcome', num2str(ni)]];
%     conditions = [conditions;['laserPrevNP', num2str(ni)]];
%     conditions =[conditions; ['laserNP', num2str(ni)]];
%     conditions = [conditions; ['laserPrevFull', num2str(ni)]];
%     conditions = [conditions; ['laserFull', num2str(ni)]];
% end
% data.conditions = [data.conditions; conditions];
data.conditions = unique(data.conditions);
%%% laser trials 
for ic = 1:numel(data.conditions)
    % Trial type indices
    clear rewIdx nrewIdx
    rewIdx = find(data.reward==1);
    nrewIdx = find(data.reward==0);
    eval(sprintf('rewIdx = rewIdx(ismember(rewIdx, data.trialIdx_%s));', data.conditions{ic}));
    eval(sprintf('nrewIdx = nrewIdx(ismember(nrewIdx, data.trialIdx_%s));', data.conditions{ic}));
    % Lever press response time
    eval(sprintf('rt.leverPress_%s = data.leverPress(data.trialIdx_%s) - data.leverPresentation(data.trialIdx_%s);', data.conditions{ic}, data.conditions{ic}, data.conditions{ic}));
    % Lever press response time by current trial type
    eval(sprintf('rt.leverRew_%s = rt.leverPress(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.leverNRew_%s = rt.leverPress(nrewIdx);', data.conditions{ic}));
    % Nose poke response times
    eval(sprintf('rt.nosePoke_%s = data.nosePokeEntry(data.trialIdx_%s)-data.trialStart(data.trialIdx_%s);', data.conditions{ic}, data.conditions{ic}, data.conditions{ic}));
    % Nose poke response times by current trial type
    eval(sprintf('rt.noseRew_%s = rt.nosePoke(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.noseNRew_%s = rt.nosePoke(nrewIdx);', data.conditions{ic}));
    % Nose port withdrawal times
    eval(sprintf('rt.withdraw_%s = data.nosePokeExit(data.trialIdx_%s) - data.nosePokeEntry(data.trialIdx_%s);', data.conditions{ic}, data.conditions{ic}, data.conditions{ic}));
    % Nose port withdrawal times by current trial type    
    eval(sprintf('rt.withdrawRew_%s = rt.withdraw(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.withdrawNRew_%s = rt.withdraw(nrewIdx);', data.conditions{ic}));
    
    
    % Previous trial type indices
    clear rewIdx nrewIdx
    rewIdx = find(data.reward==1)+1; %indices of trials following a reward
    nrewIdx = find(data.reward==0)+1; %indices of trials following no reward
    eval(sprintf('rewIdx = rewIdx(ismember(rewIdx, data.trialIdx_%s));', data.conditions{ic}));
    eval(sprintf('nrewIdx = nrewIdx(ismember(nrewIdx, data.trialIdx_%s));', data.conditions{ic})); 
    % Lever press response time by previous trial type
    eval(sprintf('rt.leverPrevRew_%s = rt.leverPress(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.leverPrevNRew_%s = rt.leverPress(nrewIdx);', data.conditions{ic}));
    % Nose poke response times by previous trial type
    eval(sprintf('rt.nosePrevRew_%s = rt.nosePoke(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.nosePrevNRew_%s = rt.nosePoke(nrewIdx);', data.conditions{ic}));
     % Nose poke withdrawal times by previous trial type
    eval(sprintf('rt.withdrawPrevRew_%s = rt.withdraw(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.withdrawPrevNRew_%s = rt.withdraw(nrewIdx);', data.conditions{ic}));
    
    clear congIdx incongIdx rewIdx nrewIdx
    congIdx = [find(data.reward==1 & data.choice == data.blocks.blockIDAll-1)+1 find(data.reward == 0 & data.choice ~= data.blocks.blockIDAll-1)+1];
    incongIdx = [find(data.reward==0 & data.choice == data.blocks.blockIDAll-1)+1 find(data.reward == 1 & data.choice ~= data.blocks.blockIDAll-1)+1];
    
    eval(sprintf('congIdx = congIdx(ismember(congIdx, data.trialIdx_%s));', data.conditions{ic}));
    eval(sprintf('incongIdx = incongIdx(ismember(incongIdx, data.trialIdx_%s));', data.conditions{ic})); 
    % Lever press response time by previous trial type
    eval(sprintf('rt.leverPrevCong_%s = rt.leverPress(congIdx);', data.conditions{ic}));
    eval(sprintf('rt.leverPrevIncong_%s = rt.leverPress(incongIdx);', data.conditions{ic}));
    % Nose poke response times by previous trial type
    eval(sprintf('rt.nosePrevCong_%s = rt.nosePoke(congIdx);', data.conditions{ic}));
    eval(sprintf('rt.nosePrevIncong_%s = rt.nosePoke(incongIdx);', data.conditions{ic}));
    
     % Nose poke withdrawal times by previous trial type
    eval(sprintf('rt.withdrawPrevCong_%s = rt.withdraw(congIdx);', data.conditions{ic}));
    eval(sprintf('rt.withdrawPrevIncong_%s = rt.withdraw(incongIdx);', data.conditions{ic}));
    

    
    % Trial duration
    rewIdx = find(data.reward==1);
    nrewIdx = find(data.reward==0);
    rewIdx = rewIdx(ismember(rewIdx+1, data.trialIdx_all));
    eval(sprintf('rewIdx = rewIdx(ismember(rewIdx,data.trialIdx_%s));', data.conditions{ic}));
    nrewIdx = nrewIdx(ismember(nrewIdx+1, data.trialIdx_all));
    eval(sprintf('nrewIdx = nrewIdx(ismember(nrewIdx,data.trialIdx_%s));', data.conditions{ic}));

    eval(sprintf('rt.trialDur_%s = data.trialStart(data.trialIdx_%s(ismember(data.trialIdx_%s+1, data.trialIdx_all))+1) - data.trialStart(data.trialIdx_%s(ismember(data.trialIdx_%s+1, data.trialIdx_all)));', data.conditions{ic},...
        data.conditions{ic}, data.conditions{ic}, data.conditions{ic}, data.conditions{ic}));
    
    eval(sprintf('rt.trialDurRew_%s = data.trialStart(rewIdx+1) - data.trialStart(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.trialDurNRew_%s = data.trialStart(nrewIdx+1) - data.trialStart(nrewIdx);', data.conditions{ic}));
    
    rewIdx = rewIdx(ismember(rewIdx+2, data.trialIdx_all))+1; 
    nrewIdx = nrewIdx(ismember(nrewIdx+2, data.trialIdx_all))+1;
    eval(sprintf('rt.trialDurPrevRew_%s = data.trialStart(rewIdx+1) - data.trialStart(rewIdx);', data.conditions{ic}));
    eval(sprintf('rt.trialDurPrevNRew_%s = data.trialStart(nrewIdx+1) - data.trialStart(nrewIdx);', data.conditions{ic}));
    
    % Reward consumption time
    clear rewIdx
    rewIdx = find(data.reward==1);
    
    eval(sprintf('rewIdx = rewIdx(ismember(rewIdx,data.trialIdx_%s));', data.conditions{ic}));
    eval(sprintf('rt.rewardConsume_%s = data.rewardExit(rewIdx) - data.rewardEntry(rewIdx);', data.conditions{ic}));
    
    
end




% Find and zscore individual sessions. . .should be separated by NaNs, catch for old files where separated by 0s

idx = find(data.choice == -10);
ii = 1;
temp = [];

if isempty(idx)
    idx = find(isnan(rt.nosePoke));
    if ~isempty(idx)
        for n = 1:2:numel(idx)
            temp = [temp zscore(rt.nosePoke(ii:idx(n)-1)) NaN NaN];
            ii = idx(n)+2;
        end
        temp(idx) = NaN;
        rt.nosePoke_zscore = temp;
    else
        rt.nosePoke_zscore = (rt.nosePoke-nanmean(rt.nosePoke))./nanstd(rt.nosePoke);
    end
elseif numel(idx)==1
    keyboard
    temp1 = rt.nosePoke;
    temp1(idx) = NaN;
    rt.nosePoke_zscore = (temp1-nanmean(temp1))./nanstd(temp1);
else
  %  keyboard
    for n = 1:2:numel(idx)
        temp = [temp zscore(rt.nosePoke(ii:idx(n)-1)) NaN NaN];
        ii = idx(n)+2;
    end
    temp(idx) = NaN;
    rt.nosePoke_zscore = temp;
end


ii = 1;
temp = [];
idx = find(rt.leverPress == 0);
if isempty(idx)
    idx = find(isnan(rt.leverPress));
    idx(data.errors(idx)~=-10) = [];
    if ~isempty(idx)
        for n = 1:2:numel(idx)
            temp = [temp (rt.leverPress(ii:idx(n)-1)-nanmean(rt.leverPress(ii:idx(n)-1)))./nanstd(rt.leverPress(ii:idx(n)-1)) NaN NaN];
            ii = idx(n)+2;
        end
       rt.leverPress_zscore = temp;      
    else
        idx = find(data.errors~=0);
        temp = rt.leverPress;
        temp(idx) = NaN;
        rt.leverPress_zscore =  (temp-nanmean(temp))./nanstd(temp);
    end
elseif numel(idx)==1
    keyboard
    temp1 = rt.leverPress;
    temp1(idx) = NaN;
    rt.leverPress_zscore =  (temp1-nanmean(temp1))./nanstd(temp1);
else
    keyboard
    ii = 1;
    temp = [];
    for n = 1:2:numel(idx)
        temp = [temp zscore(rt.leverPress(ii:idx(n)-1)) NaN NaN];
        ii = idx(n)+2;
    end
    temp(idx) = NaN;
    rt.leverPress_zscore = temp;
end

ii = 1;
temp = [];
idx = find(rt.withdraw == 0);
if isempty(idx)
    idx = find(isnan(rt.withdraw));
    if ~isempty(idx)
        for n = 1:2:numel(idx)
            temp = [temp zscore(rt.withdraw(ii:idx(n)-1)) NaN NaN];
            ii = idx(n)+2;
        end
        temp(idx) = NaN;
        rt.withdraw_zscore = temp;
    else
         rt.withdraw_zscore = (rt.withdraw-nanmean(rt.withdraw))./nanstd(rt.withdraw);
    end
elseif numel(idx)==1
    temp1 = rt.withdraw;
    temp1(idx) = NaN;
    rt.withdraw_zscore = (temp1-nanmean(temp1))./nanstd(temp1);
else
    ii = 1;
    temp = [];
    for n = 1:2:numel(idx)
        temp = [temp zscore(rt.withdraw(ii:idx(n)-1)) 0 0];
        ii = idx(n)+2;
    end
    temp(idx) = NaN;
    rt.withdraw_zscore = temp;
end

% Zscored response time for previously rewarded and unrewarded trials for data.conditions 
prevRewIdx = find(data.reward==1)+1;
prevRewIdx = prevRewIdx(ismember(prevRewIdx, data.trialIdx_all));
prevNRewIdx = find(data.reward==0)+1;
prevNRewIdx = prevNRewIdx(ismember(prevNRewIdx, data.trialIdx_all));

rt.nosePrevRew_zscore = rt.nosePoke_zscore(prevRewIdx(ismember(prevRewIdx-1, data.trialIdx)));
rt.nosePrevNRew_zscore = rt.nosePoke_zscore(prevNRewIdx(ismember(prevNRewIdx-1, data.trialIdx)));
rt.leverPrevRew_zscore = rt.leverPress_zscore(prevRewIdx(ismember(prevRewIdx-1, data.trialIdx)));
rt.leverPrevNRew_zscore = rt.leverPress_zscore(prevNRewIdx(ismember(prevNRewIdx-1, data.trialIdx)));
rt.withdrawPrevRew_zscore = rt.withdraw_zscore(prevRewIdx(ismember(prevRewIdx-1, data.trialIdx)));
rt.withdrawPrevNRew_zscore = rt.withdraw_zscore(prevNRewIdx(ismember(prevNRewIdx-1, data.trialIdx)));


for ic = 1:numel(data.conditions)
    eval(sprintf('rt.nosePrevRew_zscore_%s = rt.nosePoke_zscore(prevRewIdx(ismember(prevRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
    eval(sprintf('rt.nosePrevNRew_zscore_%s = rt.nosePoke_zscore(prevNRewIdx(ismember(prevNRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
    eval(sprintf('rt.leverPrevRew_zscore_%s = rt.leverPress_zscore(prevRewIdx(ismember(prevRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
    eval(sprintf('rt.leverPrevNRew_zscore_%s = rt.leverPress_zscore(prevNRewIdx(ismember(prevNRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
    eval(sprintf('rt.withdrawPrevRew_zscore_%s = rt.withdraw_zscore(prevRewIdx(ismember(prevRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
    eval(sprintf('rt.withdrawPrevNRew_zscore_%s = rt.withdraw_zscore(prevNRewIdx(ismember(prevNRewIdx, data.trialIdx_%s)));', data.conditions{ic}, data.conditions{ic}));
end

% Find nose poke and lever press times following multiple trials of a given type in a row

for ri = 1:length(data.multiRewIdx)
    data.multiRewIdx{ri} = data.multiRewIdx{ri}(ismember(data.multiRewIdx{ri}+1, data.trialIdx_all));
    rt.multiPrevRewNose{ri} = rt.nosePoke(data.multiRewIdx{ri}+1);
    rt.multiPrevRewLever{ri} = rt.leverPress(data.multiRewIdx{ri}+1);
end

for ri = 1:length(data.multiNrewIdx)
    data.multiNrewIdx{ri} = data.multiNrewIdx{ri}(ismember(data.multiNrewIdx{ri}+1, data.trialIdx_all));
    rt.multiPrevNrewNose{ri} = rt.nosePoke(data.multiNrewIdx{ri}+1);
    rt.multiPrevNrewLever{ri} = rt.leverPress(data.multiNrewIdx{ri}+1);
end



% For multiple trials in a row (zscore)
for ri = 1:length(data.multiRewIdx)
    rt.multiPrevRewNose_zscore{ri} = rt.nosePoke_zscore(data.multiRewIdx{ri}+1);
    rt.multiPrevRewLever_zscore{ri} = rt.leverPress_zscore(data.multiRewIdx{ri}+1);
end

for ri = 1:length(data.multiNrewIdx)
    rt.multiPrevNrewNose_zscore{ri} = rt.nosePoke_zscore(data.multiNrewIdx{ri}+1);
    rt.multiPrevNrewLever_zscore{ri} = rt.leverPress_zscore(data.multiNrewIdx{ri}+1);
end

% Laser trials 
for ri = 1:length(data.multiRewIdx)
    for ne = 1:numel(data.conditions)
        clear temp
        temp = eval(sprintf('data.multiRewIdx{ri}(ismember(data.multiRewIdx{ri}+1, data.trialIdx_%s));', data.conditions{ne}));
        eval(sprintf('rt.multiPrevRewNose_%s{ri} = rt.nosePoke(temp+1);', data.conditions{ne}));
        eval(sprintf('rt.multiPrevRewLever_%s{ri} = rt.leverPress(temp+1);', data.conditions{ne}));
        eval(sprintf('rt.multiPrevRewNose_zscore_%s{ri} = rt.nosePoke_zscore(temp+1);', data.conditions{ne}));
        eval(sprintf('rt.multiPrevRewLever_zscore_%s{ri} = rt.leverPress_zscore(temp+1);', data.conditions{ne}));
    end
end


for ri = 1:length(data.multiNrewIdx)
     for ne = 1:numel(data.conditions)
    temp = eval(sprintf('data.multiNrewIdx{ri}(ismember(data.multiNrewIdx{ri}+1, data.trialIdx_%s));', data.conditions{ne}));
    eval(sprintf('rt.multiPrevNRewNose_%s{ri} = rt.nosePoke(temp+1);', data.conditions{ne}));
    eval(sprintf('rt.multiPrevNrewLever_%s{ri} = rt.leverPress(temp+1);', data.conditions{ne}));
    eval(sprintf('rt.multiPrevNrewNose_zscore_%s{ri} = rt.nosePoke_zscore(temp+1);', data.conditions{ne}));
    eval(sprintf('rt.multiPrevNrewLever_zscore_%s{ri} = rt.leverPress_zscore(temp+1);', data.conditions{ne}));
      end
end

for ri = 1:length(data.multiRewIdx)
    temp = data.multiRewIdx{ri}(ismember(data.multiRewIdx{ri}, data.trialIdx));
    temp = temp(ismember(temp+1, data.trialIdx));  
    rt.multiPrevRewNose_noLaser{ri} = rt.nosePoke(temp+1);
    rt.multiPrevRewLever_noLaser{ri} = rt.leverPress(temp+1);
    rt.multiPrevRewNose_zscore_noLaser{ri} = rt.nosePoke_zscore(temp+1);
    rt.multiPrevRewLever_zscore_noLaser{ri} = rt.leverPress_zscore(temp+1);
end

for ri = 1:length(data.multiNrewIdx)
    temp = data.multiNrewIdx{ri}(ismember(data.multiNrewIdx{ri}, data.trialIdx));
    temp = temp(ismember(temp+1, data.trialIdx));
    rt.multiPrevNrewNose_noLaser{ri} = rt.nosePoke(temp+1);
    rt.multiPrevNrewLever_noLaser{ri} = rt.leverPress(temp+1);
    rt.multiPrevNrewNose_zscore_noLaser{ri} = rt.nosePoke_zscore(temp+1);
    rt.multiPrevNrewLever_zscore_noLaser{ri} = rt.leverPress_zscore(temp+1);
end

%%% Find nose poke and lever press times following multiple laser trials in a row
conditions = data.conditions;
idx = strcmp(conditions, 'all');
conditions = conditions(~idx);

for ic = 1:numel(conditions)
     if isempty(eval(sprintf('(data.trialIdx_multi%s)', conditions{ic})))
          eval(sprintf('rt.noseMulti%s{1} = [];', conditions{ic}));
        eval(sprintf('rt.leverMulti%s{1} = [];', conditions{ic}));
        eval(sprintf('rt.noseMulti%s_zscore{1} =[];', conditions{ic}));
        eval(sprintf('rt.leverMulti%s_zscore{1} =[];', conditions{ic}));
     else
    for ri = 1:eval(sprintf('length(data.trialIdx_multi%s)', conditions{ic}))
        eval(sprintf('rt.noseMulti%s{ri} = rt.nosePoke(data.trialIdx_multi%s{ri});', conditions{ic}, conditions{ic}));
        eval(sprintf('rt.leverMulti%s{ri} = rt.leverPress(data.trialIdx_multi%s{ri});', conditions{ic}, conditions{ic}));
        eval(sprintf('rt.noseMulti%s_zscore{ri} = rt.nosePoke_zscore(data.trialIdx_multi%s{ri});', conditions{ic}, conditions{ic}));
        eval(sprintf('rt.leverMulti%s_zscore{ri} = rt.leverPress_zscore(data.trialIdx_multi%s{ri});', conditions{ic}, conditions{ic}));
    end
     end
end




data.rt = rt;
