function dataForQ_humanBandit(cohorts,cohort)

basefilename = fullfile(whereAreWe('figurecode'),'raw_data','figure_1','human_bandit'); % where is the data
% cohorts      = {'cohort1';'cohort2'}; % which data to use
% cohort       = 'all'; %for saving

%%% Configure data to run Qlearning model in stan 

sList = [];
for nc = 1:numel(cohorts)
    temp = dir(fullfile(basefilename, [cohorts{nc} '_f'], ['*.mat']));
    for ns = 1:numel(temp)
        sList = cat(1,sList,{fullfile(temp(ns).folder, temp(ns).name)}); 
    end
    temp = dir(fullfile(basefilename, [cohorts{nc} '_m'], ['*.mat']));
    for ns = 1:numel(temp)
        sList = cat(1,sList,{fullfile(temp(ns).folder, temp(ns).name)}); 
    end
end

%trial number, subject number, action (right = 0, left = 1), reward, session, group, (trial start?)

Trial = [];
SubjID = [];
Action = [];
Reward = [];
exclude = [];

trialStartLatency = [];
choiceLatency = []; 

sCounter = 1;

    
for na = 1:numel(sList)
    load(sList{na})
    stay = [NaN; data.choice(1:end-1)==data.choice(2:end)];
    prevReward = [NaN; data.reward(1:end-1)];
    if sum(~isnan(data.reward))<100
        exclude = cat(1,exclude,na); 
 
    else
        
    Trial  = cat(1,Trial,[1:sum(~isnan(data.reward))]');
    Reward = cat(1, Reward, data.reward(~isnan(data.reward)));
    SubjID = cat(1, SubjID, ones(sum(~isnan(data.reward)),1).*sCounter);
    Action = cat(1, Action, data.choice(~isnan(data.reward)));
    
    % latencies
    if ~isempty(data.trialStartLatency)
    trialStartLatency = cat(1,trialStartLatency,data.trialStartLatency(~isnan(data.reward)));
    else
        trialStartLatency = cat(1, trialStartLatency, nan(size(data.reward(~isnan(data.reward)))));
    end
    choiceLatency     = cat(1,choiceLatency,data.choiceLatency(~isnan(data.reward))); 
    
    sCounter = sCounter+1;
    fprintf('%s \n',num2str(sCounter-1))
    end
end

sList = sList(~ismember(1:numel(sList),exclude));
NT = max(Trial); % max number of trials (per subject)
NS = numel(sList); % number of subjects
NT_all = zeros(NS,1); % number of trials per subject
for n = 1:numel(sList)
    NT_all(n) = max(Trial(SubjID==n));
end


r = ones(NS,NT).*-1;
c = ones(NS,NT).*-1; 

for n = 1:numel(sList)
        c(n,1:sum(SubjID==n)) = Action(SubjID==n);
        r(n,1:sum(SubjID==n)) = Reward(SubjID==n);
end
c(isnan(r)) = -1;
r(isnan(r)) = -1;

    
data = [SubjID,Trial,Action,Reward];
savehere = fullfile(whereAreWe('figurecode'),'processed_data','dataForQ','human_bandit',cohort);
mkdir(savehere); 
writematrix(data,fullfile(savehere,'dataForQ.csv'));
save(fullfile(savehere,'sList'), 'sList')
save(fullfile(savehere,'data_toStan.mat'),'NT','NS','NT_all','c','r','trialStartLatency','choiceLatency')

