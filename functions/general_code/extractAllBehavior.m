function [behaviorTable, daylist] = extractAllBehavior(id,perfThresh,qFile,sheetnames,weightSpreadsheet,ext_laser,femaleFlag,binNum,binNum_choice,cohort,sessionLength,intervalThresh) 
 
% Extract session by session behavior data

%%% Inputs %%%
%id: animal ID (string)
%perfThresh: performance threshold for inclusion 
%qFile: which Q-learning model 
%sheetnames: names of the pages in the weight spreadsheet 
%weightSpreadsheet: location of the weight spreadsheet
%ext_laser: file extension for laser sessions 
%femaleFlag: is this a female? (0 or 1)

%%% Outputs %%%
%behaviorTable.choice: which lever (right=0,left=1)
%behaviorTable.reward: reward=1, no reward = 0
%behaviorTable.previousReward: reward=1,no reward=0
%behaviorTable.trialInit: trial initiation latency (seconds)
%behaviorTable.trialInit_zscore: session by session zscored trial initiation latency 
%behaviorTable.leverPressLatency: lever press latency (seconds) 
%behaviorTable.leverPressLatency_zscore: session by session zscored lever press latency 
%behaviorTable.session = session indicator (not sorted) 
%behaviorTable.trial = trial number
%behaviorTable.datenum = date of session
%behaviorTable.weight = animal's post-training weight 
%behaviorTable.values = various q-values 
%behaviorTable.female = is the animal female? 
%behaviorTable.laser  = laser trial? 1 = previous outcome, 2 = previous nose poke, 3 = ITI 
%behaviorTable.aID    = animal ID; 
%behaviorTable.laserSession = is this a laser session? 
%% Initialize variables for output table
choice            = [];
reward            = [];
previousReward    = [];
trialInit         = [];
trialInit_zscore  = [];
trialInit_thresh  = []; 
leverPress        = [];
leverPress_zscore = [];
session           = [];
trial             = [];
sessionDate       = [];
weight            = [];
female            = [];
laser             = []; 
laserSession      = [];
aID               = []; 
stay              = []; 
trialTime         = [];
block             = [];
ipsiChoice        = []; 
valTypes = {'qTot';'qChosenDiff';'qDiff'};

for nv = 1:numel(valTypes)
    eval(sprintf('%s = [];',valTypes{nv}))
end

%% Parameters 

basefilename  = fullfile(whereAreWe('behavior'), id);
% generate list of animal IDs from opto experiment 
if ~contains(cohort,'d1') && ~contains(cohort,'d2')
ids_cohort = cat(1,generateAnimalList(sprintf('%s_female',cohort)), generateAnimalList(sprintf('%s_male',cohort)),generateAnimalList(sprintf('%s_yfp_female',cohort)), generateAnimalList(sprintf('%s_yfp_male',cohort))); 
else
    ids_cohort = cat(1,generateAnimalList(sprintf('%s_female',cohort)), generateAnimalList(sprintf('%s_male',cohort)),generateAnimalList('DMS_yfp_female'), generateAnimalList('DMS_yfp_male')); 
end
%% Find weights for animal id on spreadsheet and extract dates and values 
sheetIdx = find(contains(sheetnames,['_' id '_']));
if isempty(sheetIdx)
    sheetIdx = find(contains(sheetnames, ['_' id]));
end
try
    opts = detectImportOptions(weightSpreadsheet,'Sheet',sheetIdx);
    aIdx = find(contains(opts.SelectedVariableNames,id));
    weightData = readmatrix(weightSpreadsheet, 'Sheet',(sheetIdx),'OutputType','string');
    if isempty(aIdx)
        aIdx = find(contains(weightData(1,:),id));
        weightData = weightData(2:end,:);
    end
    weights = double(weightData(:,aIdx));
    weightDates = cell2mat(arrayfun(@(x) datenum(x),weightData(:,aIdx-1),'UniformOutput',false));
    temp = nan(size(weights));
    try
        temp(~isnan(weights)) = weightDates;
    catch
        temp = weightDates;
    end
    weightDates = temp;
    cageWeights = double(weightData(:,2:2:end));
catch
    keyboard
    warning(sprintf('No weights for %s',id))
end
    
%% Load session list and values for control data
keyboard
try dataAll = load(fullfile(basefilename, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)));
catch
    valueExtraction_TAB({id},qFile,sessionLength,perfThresh);
    dataAll = load(fullfile(basefilename, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)));
end
%% Load data for laser sessions 

if sum(strcmp(ids_cohort,id))>0
    try dataOpto = load(fullfile(basefilename, sprintf('value%s_perfThresh%s_%s',ext_laser,num2str(0),'qLearn_session_all.mat')));
    catch
        valueExtraction_opto({id},qFile,ext_laser,0);
        dataOpto = load(fullfile(basefilename, sprintf('value%s_perfThresh%s_%s',ext_laser,num2str(0),'qLearn_session_all.mat')));
  end

%% concatenate laser data 
sessionIdx = unique(dataOpto.session);
sessionCounter = 1;
stay              = cat(1, stay, dataOpto.stay');
choice            = cat(1,choice,dataOpto.leftChoice');
ipsiChoice        = cat(1,ipsiChoice,dataOpto.leftChoice'==dataOpto.laserSide);
previousReward    = cat(1,previousReward,dataOpto.prevOutcome');
trialInit         = cat(1,trialInit,dataOpto.trialStart');
trialInit_zscore  = cat(1,trialInit_zscore,dataOpto.trialStart_zscore');
temp              = nan(size(dataOpto.trialStart'));
tempSession       = unique(dataOpto.session);
for ns = 1:numel(tempSession)
    thisIdx = find(dataOpto.session==tempSession(ns));
    thisLatency = dataOpto.trialStart(thisIdx);
    stopIdx = find(thisLatency>intervalThresh,1,'first');
    if isempty(stopIdx)
        temp(thisIdx) = thisLatency;
    else
        temp(thisIdx(1:stopIdx-1)) = thisLatency(1:stopIdx-1);
    end
end
trialInit_thresh  = cat(1,trialInit_thresh,temp); 
leverPress        = cat(1,leverPress,dataOpto.leverPress');
leverPress_zscore = cat(1,leverPress_zscore,dataOpto.leverPress_zscore');
female            = cat(1,female,ones(size(dataOpto.leftChoice')).*femaleFlag);


laserType         = dataOpto.laserType;
laserType(laserType==-1) = 0;
laserType(laserType==6) = 3;
laserIdx1          = find(laserType==1);
laserIdx2          = find(laserType==2);
laserIdx3          = find(laserType==3);
laserType(laserIdx1)= -1; % trial w/ laser
laserType(laserIdx2)= -2; % trial w/ laser
laserType(laserIdx2+1) = 2; % trial after laser
laserType(laserIdx1+1) = 1; % trial after laser
laserType(laserIdx3+1)= -3; % trial after laser (for ITI condition)
laserType(laserIdx3) = 3; 
laser             = cat(1,laser,laserType(1:numel(dataOpto.stay))');

laserSession      = cat(1,laserSession,ones(size(dataOpto.leftChoice')));
aID               = cat(1,aID,repmat({id},size(dataOpto.leftChoice')));


for nv = 1:numel(valTypes)
    eval(sprintf('%s = cat(1,%s,dataOpto.%s);',valTypes{nv},valTypes{nv},valTypes{nv}))
end

for nf = 1:numel(dataOpto.flist)
    %try
        load(fullfile(basefilename,[dataOpto.flist{nf}(2:end-4),sprintf('_%s.mat',ext_laser)]));
        data.trialIdx_all = 1:numel(data.choice);
        reward        = cat(1,reward,data.reward(data.trialIdx_all)',-10);
        sessionDate   = cat(1,sessionDate,(ones(size(data.trialIdx_all'))).*datenum(data.raw.header.Start_Date),-10);
        block = cat(1, block,data.blocks.blockIDAll(data.trialIdx_all)',-10);
        try
            weight        = cat(1,weight,ones(size(data.trialIdx_all')).*weights(weightDates==datenum(data.raw.header.Start_Date)),-10);
        catch
            weight        = cat(1,weight,nan(size(data.trialIdx_all')),-10);
        end
        trial         = cat(1,trial,data.trialIdx_all',-10);
        session       = cat(1,session,ones(size(data.trialIdx_all')).*sessionCounter,-10);
        trialTime     = cat(1,trialTime, data.trialStart(data.trialIdx_all)',-10);
        
        sessionCounter= sessionCounter+1;
%     catch
%         keyboard
%         reward        = cat(1,reward,nan(size(data.reward(data.trialIdx_all)')));
%         sessionDate   = cat(1,sessionDate,nan(size(data.trialIdx_all')));
%         try
%             weight        = cat(1,weight,nan(size(data.trialIdx_all')));
%         catch
%             weight        = cat(1,weight,nan(size(data.trialIdx_all')));
%         end
%         trial         = cat(1,trial,nan(size(data.trialIdx_all')));
%         session       = cat(1,session,nan(size(data.trialIdx_all')));
%         %sessionCounter= sessionCounter+1;
%     end
    
    
end
%%
%     
%     load(fullfile(basefilename, sprintf('laserSummary_%s.mat',ext_laser)));
%     if ~isfield(data,'qFile')
%         Load Q-values and extract laser sessions
%         load(fullfile(basefilename, qFile));
%         data = concatValue(data,qLearn,id);
%         data.qFile = qFile;
%         save(fullfile(basefilename, ['laserSummary_' ext_laser '.mat']), 'data')
%    elseif ~contains(data.qFile,qFile)
%         Load Q-values and extract laser sessions
%         load(fullfile(basefilename,  qFile));
%         data = concatValue(data,qLearn,id);
%         data.qFile = qFile;
%         save(fullfile(basefilename,  ['laserSummary_' ext_laser '.mat']), 'data')
%     end
%     
%     
%     % Standardize filename format for laser and non-laser sessions
%     idx     = cellfun(@(x) strfind(x,'/'), data.daylist,'UniformOutput',false);
%     if isempty(idx{1})
%         idx     = cellfun(@(x) strfind(x,'\'), data.daylist,'UniformOutput',false);
%     end
%     idx         = cellfun(@(x) x(end)+1, idx,'UniformOutput',false);
%     daylist     = cellfun(@(x,y) x(y:end),data.daylist,idx,'UniformOutput',false);
%     tempDir     = dir(fullfile(basefilename, sprintf('*.%s',ext_laser)));
%     tempDir     = {tempDir(:).name};
%     daylist_new = cellfun(@(x) tempDir(contains(tempDir,x(1:end-8))),daylist,'UniformOutput',false);
%     
%     daylist = cat(1,daylist_new',dataAll.flist);
%     if size(daylist_new,1) == 1 && size(daylist_new{1},1)==0
%         
%         daylist           = dataAll.flist;
%         daylist_laser     = [];
%         sessionCounter    = 1;
%     else
%         % Concatenate laser data
%         sessionBreaks = find(data.choice==-10);
%         stay              = cat(1,stay,data.choice(data.trialIdx_all)'==data.choice(data.trialIdx_all-1)');
%         choice            = cat(1,choice,data.choice(data.trialIdx_all)');
%         reward            = cat(1,reward,data.reward(data.trialIdx_all)');
%         prevReward_temp   = cat(1,NaN,data.reward(1:end-1)');
%         previousReward    = cat(1,previousReward,prevReward_temp(data.trialIdx_all));
%         trialInit         = cat(1,trialInit,data.rt.nosePoke(data.trialIdx_all)');
%         temp              = cell2mat(arrayfun(@(x,y) [nanzscore(data.rt.nosePoke(x:y)) -10 -10], [1 sessionBreaks(2:2:end-1)+1], [sessionBreaks(1:2:end)-1],'UniformOutput',false));
%         trialInit_zscore  = cat(1,trialInit_zscore,temp(data.trialIdx_all)');
%         leverPress        = cat(1,leverPress,data.rt.leverPress(data.trialIdx_all)');
%         temp              = cell2mat(arrayfun(@(x,y) [nanzscore(data.rt.leverPress(x:y)) -10 -10], [1 sessionBreaks(2:2:end)+1], [sessionBreaks(1:2:end)-1 numel(data.choice)],'UniformOutput',false));
%         leverPress_zscore = cat(1,leverPress_zscore,temp(data.trialIdx_all)');
%         temp              = cell2mat(arrayfun(@(x,y,z) [ones(size(x:y)).*z -10 -10], [1 sessionBreaks(2:2:end-1)+1], [sessionBreaks(1:2:end)-1],1:numel(data.daylist),'UniformOutput',false));
%         session           = cat(1,session,temp(data.trialIdx_all)');
%         temp              = cell2mat(arrayfun(@(x,y) [1:numel(x:y) -10 -10], [1 sessionBreaks(2:2:end-1)+1], [sessionBreaks(1:2:end)-1],'UniformOutput',false));
%         trial             = cat(1,trial,temp(data.trialIdx_all)');
%         female            = cat(1,female,ones(size(data.trialIdx_all))'.*femaleFlag);
%         aID               = cat(1,repmat({id},size(data.trialIdx_all')));
%         % laser types (1=outcome, 2= nose poke, 3=ITI)
%         laserType         = data.laserType;
%         laserType(laserType==-1) = 0;
%         laserType(laserType==6) = 3;
%         laserIdx          = find(laserType==1);
%         laserType(laserIdx)= -1; % trial w/ laser
%         laserType(laserIdx+1) = 1; % trial after laser
%         laserIdx          = find(laserType==2);
%         laserType(laserIdx)= -2; % trial w/ laser
%         laserType(laserIdx+1) = 2; % trial after laser
%         laserIdx          = find(laserType==3);
%         laserType(laserIdx+1)= -3; % trial after laser (for ITI condition)
%         laser             = cat(1,laser,laserType(data.trialIdx_all)');
%         laserSession      = cat(1,laserSession,ones(size(data.trialIdx_all))');
%         % concatenate values
%         data.qTot = data.QRight'+data.QLeft';
%         data.qDiff = data.QRight'-data.QLeft';
%         data.qChosenDiff = nan(size(data.QRight))';
%         data.qChosenDiff(data.choice==0) = data.QRight(data.choice==0)-data.QLeft(data.choice==0);
%         data.qChosenDiff(data.choice==1) = data.QLeft(data.choice==1)-data.QRight(data.choice==1);
%         
%         for nv = 1:numel(valTypes)
%             try
%                 eval(sprintf('%s = cat(2,%s,data.%s(data.trialIdx_all));',valTypes{nv},valTypes{nv},valTypes{nv}));
%             catch
%                 keyboard
%             end
%         end
%         
%         % load each data file to extract session date and find session weight
%         daylist_laser = data.daylist;
%         tempDate   = [];
%         tempWeight = [];
%         for nf = 1:numel(daylist_laser)
%             idx = strfind(daylist_laser{nf},'\');
%             if isempty(idx)
%                 idx = strfind(daylist_laser{nf}, '/');
%             end
%             load(fullfile(basefilename,daylist_laser{nf}(idx(end)+1:end)));
%             tempDate = cat(1,tempDate,(ones(size(data.trialIdx_all'))).*datenum(data.raw.header.Start_Date));
%             try
%                 tempWeight = cat(1,tempWeight,ones(size(data.trialIdx_all')).*weights(weightDates==datenum(data.raw.header.Start_Date)));
%             catch
%                 tempWeight        = cat(1,tempWeight,nan(size(data.trialIdx_all')));
%             end
%         end
%         
%         sessionDate       = cat(1,sessionDate,tempDate);
%         weight            = cat(1,weight,tempWeight);
%         
%         
%         sessionCounter    = numel(daylist_laser)+1;
%     end
%     clear daylist_new tempDir idx
%    
daylist = dataOpto.flist;
else
    % if no laser sessions
    daylist           = dataAll.flist;
    daylist_laser     = [];
    sessionCounter    = 1;
end


%% concatenate data for non laser sessions 
sessionIdx = unique(dataAll.session);

stay              = cat(1, stay, dataAll.stay');
choice            = cat(1,choice,dataAll.leftChoice');
ipsiChoice        = cat(1,ipsiChoice,dataAll.leftChoice'); 
previousReward    = cat(1,previousReward,dataAll.prevOutcome');
trialInit         = cat(1,trialInit,dataAll.trialStart');
trialInit_zscore  = cat(1,trialInit_zscore,dataAll.trialStart_zscore');
temp              = nan(size(dataAll.trialStart'));
tempSession       = unique(dataAll.session);
for ns = 1:numel(tempSession)
    thisIdx = find(dataAll.session==tempSession(ns));
    thisLatency = dataAll.trialStart(thisIdx);
    stopIdx = find(thisLatency>intervalThresh,1,'first');
    if isempty(stopIdx)
        temp(thisIdx) = thisLatency;
    else
        temp(thisIdx(1:stopIdx-1)) = thisLatency(1:stopIdx-1);
    end
end
trialInit_thresh  = cat(1,trialInit_thresh,temp); 
leverPress        = cat(1,leverPress,dataAll.leverPress');
leverPress_zscore = cat(1,leverPress_zscore,dataAll.leverPress_zscore');
female            = cat(1,female,ones(size(dataAll.leftChoice')).*femaleFlag);
laser             = cat(1,laser,zeros(size(dataAll.leftChoice')));
laserSession      = cat(1,laserSession,zeros(size(dataAll.leftChoice')));
aID               = cat(1,aID,repmat({id},size(dataAll.leftChoice')));


valTypes = {'qTot';'qChosenDiff';'qDiff'};

for nv = 1:numel(valTypes)
    eval(sprintf('%s = cat(1,%s,dataAll.%s);',valTypes{nv},valTypes{nv},valTypes{nv}))
end

for nf = 1:numel(dataAll.flist)
    try
        load(fullfile(basefilename,[dataAll.flist{nf}(2:end-4),'_TAB.mat']));
        reward        = cat(1,reward,data.reward(data.trialIdx_all)');
        sessionDate   = cat(1,sessionDate,(ones(size(data.trialIdx_all'))).*datenum(data.raw.header.Start_Date));
        try
            weight        = cat(1,weight,ones(size(data.trialIdx_all')).*weights(weightDates==datenum(data.raw.header.Start_Date)));
        catch
            weight        = cat(1,weight,nan(size(data.trialIdx_all')));
        end
        trial         = cat(1,trial,data.trialIdx_all');
        session       = cat(1,session,ones(size(data.trialIdx_all')).*sessionCounter);
        trialTime     = cat(1,trialTime, data.trialStart(data.trialIdx_all)');
        block         = cat(1, block,data.blocks.blockIDAll(data.trialIdx_all)');
        sessionCounter= sessionCounter+1;
    catch
        keyboard
        reward        = cat(1,reward,nan(size(data.reward(data.trialIdx_all)')));
        sessionDate   = cat(1,sessionDate,nan(size(data.trialIdx_all')));
        try
            weight        = cat(1,weight,nan(size(data.trialIdx_all')));
        catch
            weight        = cat(1,weight,nan(size(data.trialIdx_all')));
        end
        trial         = cat(1,trial,nan(size(data.trialIdx_all')));
        session       = cat(1,session,nan(size(data.trialIdx_all')));
        trialTime     = cat(1,trialTime, data.trialStart(data.trialIdx_all)');
        block         = cat(1,block,nan(size(data.reward(data.trialIdx_all)')));


        %sessionCounter= sessionCounter+1;
    end
    
    
end
daylist = cat(1,daylist,dataAll.flist); 

%% calculate percentiles for 

for nv = 1:numel(valTypes)
    ptiles = prctile(eval(sprintf('%s',valTypes{nv})),linspace(0,100,binNum+1));
    try
    [~,~,value_ptile] = histcounts(eval(sprintf('%s',valTypes{nv})),ptiles);
    catch
        keyboard
    end
    eval(sprintf('%s_quant = value_ptile;',valTypes{nv}));
    ptiles = prctile(eval(sprintf('%s',valTypes{nv})),linspace(0,100,binNum_choice+1));
    [~,~,value_ptile] = histcounts(eval(sprintf('%s',valTypes{nv})),ptiles);
    eval(sprintf('%s_quant_choice = value_ptile;',valTypes{nv}));
    % percentile only until threshold 
    ptiles = prctile(eval(sprintf('%s(~isnan(trialInit_thresh))',valTypes{nv})),linspace(0,100,binNum+1));
    try
    [~,~,value_ptile] = histcounts(eval(sprintf('%s',valTypes{nv})),ptiles);
    catch
        keyboard
    end
    eval(sprintf('%s_quant_thresh = value_ptile;',valTypes{nv}));
    ptiles = prctile(eval(sprintf('%s(~isnan(trialInit_thresh))',valTypes{nv})),linspace(0,100,binNum_choice+1));
    [~,~,value_ptile] = histcounts(eval(sprintf('%s',valTypes{nv})),ptiles);
    eval(sprintf('%s_quant_choice_thresh = value_ptile;',valTypes{nv}));
end

try
behaviorTable                   = table(choice,'VariableNames',{'choice'});
behaviorTable.reward            = reward;
behaviorTable.previousReward    = previousReward;
behaviorTable.trialInit         = trialInit;
behaviorTable.trialInit_zscore  = trialInit_zscore;
behaviorTable.trialInit_thresh  = trialInit_thresh; 
behaviorTable.leverPress        = leverPress;
behaviorTable.leverPress_zscore = leverPress_zscore;
behaviorTable.session           = session;
behaviorTable.trial             = trial;
behaviorTable.sessionDate       = sessionDate;
behaviorTable.weight            = weight;
behaviorTable.female            = female;
behaviorTable.laser             = laser; 
behaviorTable.laserSession      = laserSession;
behaviorTable.aID               = aID; 
behaviorTable.stay              = stay;
behaviorTable.trialTime         = trialTime;
behaviorTable.block             = block;
behaviorTable.ipsiChoice        = ipsiChoice;
for nv = 1:numel(valTypes)
   eval(sprintf('behaviorTable.%s_quant = %s_quant;',valTypes{nv},valTypes{nv})); 
   eval(sprintf('behaviorTable.%s_quant_choice = %s_quant_choice;',valTypes{nv},valTypes{nv}));
   eval(sprintf('behaviorTable.%s_thresh_quant = %s_quant_thresh;',valTypes{nv},valTypes{nv})); 
   eval(sprintf('behaviorTable.%s_thresh_quant_choice = %s_quant_choice_thresh;',valTypes{nv},valTypes{nv}));
   eval(sprintf('behaviorTable.%s = %s;',valTypes{nv},valTypes{nv})); 
end
behaviorTable(behaviorTable.choice==-10,:) = [];
behaviorTable(behaviorTable.choice==-1,:) = [];
catch
    keyboard
end




end






