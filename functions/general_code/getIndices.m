function data = getIndices(data) 
%extract trial type indices for "data" .mat files (i.e. not "summary") 

%%% Get previously rewarded and unrewarded trial indices for all trial types in data.conditions

%%% non-laser trials %%%
% previously rewarded non-laser trials (both rewarded trial and following reward are non-laser trials)
idx.prevRew = find(data.reward==1);
idx.prevRew = idx.prevRew(ismember(idx.prevRew, data.trialIdx))+1;
idx.prevRew = idx.prevRew(ismember(idx.prevRew, data.trialIdx));
% previously unrewarded non-laser trials (both rewarded trial and following reward are non-laser trials)
idx.prevNRew = find(data.reward==0 & data.choice~=-1);
idx.prevNRew = idx.prevNRew(ismember(idx.prevNRew, data.trialIdx))+1;
idx.prevNRew = idx.prevNRew(ismember(idx.prevNRew, data.trialIdx));

for ne = 1:numel(data.conditions)
   eval(sprintf('idx.prevRew_%s = find(data.reward==1)+1;', data.conditions{ne}));
   eval(sprintf('idx.prevRew_%s = idx.prevRew_%s(ismember(idx.prevRew_%s, data.trialIdx_%s));', data.conditions{ne},data.conditions{ne},data.conditions{ne},data.conditions{ne}));  
   eval(sprintf('idx.prevNRew_%s = find(data.reward==0 & data.choice~=-1)+1;', data.conditions{ne}));
   eval(sprintf('idx.prevNRew_%s = idx.prevNRew_%s(ismember(idx.prevNRew_%s, data.trialIdx_%s));', data.conditions{ne},data.conditions{ne},data.conditions{ne},data.conditions{ne}));  
end

%%% Get previously ipsi and contra trial indices for all trial types in data.conditions

%%% non-laser trials %%%
% previously ipsi non-laser trials
idx.prevIpsi = find(data.ipsiChoice==1);
idx.prevIpsi = idx.prevIpsi(ismember(idx.prevIpsi, data.trialIdx))+1;
idx.prevIpsi = idx.prevIpsi(ismember(idx.prevIpsi, data.trialIdx));
% previously contra non-laser trials
idx.prevContra = find(data.ipsiChoice==0);
idx.prevContra = idx.prevContra(ismember(idx.prevContra, data.trialIdx))+1;
idx.prevContra = idx.prevContra(ismember(idx.prevContra, data.trialIdx));

for ne = 1:numel(data.conditions)
   eval(sprintf('idx.prevIpsi_%s = find(data.ipsiChoice==1)+1;', data.conditions{ne}));
   eval(sprintf('idx.prevIpsi_%s = idx.prevIpsi_%s(ismember(idx.prevIpsi_%s, data.trialIdx_%s));', data.conditions{ne},data.conditions{ne},data.conditions{ne},data.conditions{ne}));  
   eval(sprintf('idx.prevContra_%s = find(data.ipsiChoice==0)+1;', data.conditions{ne}));
   eval(sprintf('idx.prevContra_%s = idx.prevContra_%s(ismember(idx.prevContra_%s, data.trialIdx_%s));', data.conditions{ne},data.conditions{ne},data.conditions{ne},data.conditions{ne}));  
end



data.idx = idx;