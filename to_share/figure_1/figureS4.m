function stats=figureS4(extractFlag,qFile,sessionLength,perfThresh,binNum,ids)

%% Parameters
if nargin == 0
    extractFlag     = 0; 
    sessionLength   = 'long';
    perfThresh      = .1;
    qFile           = 'qLearn_session_estrous_2022.mat';
    binNum          = 4;
    binNum_choice   = 9;
    ids = generateAnimalList('estrous');
elseif nargin < 2
    sessionLength   = 'long';
    perfThresh      = .1;
    qFile           = 'qLearn_session_estrous_2022.mat';
    binNum          = 4;
    ids = generateAnimalList('estrous');
elseif nargin < 6
    ids = generateAnimalList('estrous');
end

intervalThresh =300;
basefilename = whereAreWe('behavior');
sessionStart = datenum('3/1/22');
stopDate = datenum('5/15/2022');
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

[~,fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC);
fmap = fmap(round(linspace(20,size(fmap,1),3)),:);
%% import estrous data

% % google doc info 
% ID = '1vkJ-LMZPUq7_1Fv-Y-ksp-4trC_08-0QIWL4IrRFSx8';
% 
% sheet_name = 'weight';
% url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
%     ID, sheet_name);
% sheet_name = 'estrous_consensus';
% url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
%     ID, sheet_name);
% sheet_name = 'estrous_strict';
% url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
%     ID, sheet_name);

weight_data = readtable(fullfile(whereAreWe('figurecode'),'raw_data','weight.csv'));

estrous_consensus = readtable(fullfile(whereAreWe('figurecode'),'raw_data','estrous_consensus.csv'));

estrous_strict = readtable(fullfile(whereAreWe('figurecode'),'raw_data','estrous_strict.csv'));

%% Extract q-values and latencies for each animal 
if extractFlag || ~exist(fullfile(whereAreWe('figurecode'),'processed_data','estrous',sprintf('estrousData_%s_%s',num2str(perfThresh),qFile)),'file')
    Animal = [];
    Latency = [];
    Choice = [];
    Value = [];
    Value_choice = [];
    Value_ptile = [];
    Value_ptile_choice = [];
    Estrous = [];
    Estrous_strict = [];
    Weight = [];
    Reward = [];
    Session = [];
    Latency_zscore = [];
    TrialTime = [];
    Dates = [];
    Trial = []; 
    
    dates = estrous_consensus(:,1);
    dates = datenum(dates.Date);
    dates_strict = estrous_strict(:,1);
    dates_strict = datenum(dates_strict.Date);
    weight_data = weight_data(3:end,:);
    dates_weight = weight_data(:,1);
    dates_weight = datenum(dates_weight.Date);
    for na = 1:numel(ids)
        
        try load(fullfile(basefilename,ids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
            load(fullfile(basefilename,ids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
            thisFlist{na} = flist;
        catch
            valueExtraction_TAB(ids(na),qFile,sessionLength,perfThresh);
            load(fullfile(basefilename,ids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
            load(fullfile(basefilename,ids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
            thisFlist{na} = flist;
        end
        
        
        %initialize variables
        thisLatency = [];
        thisValue   = [];
        thisValue_choice = [];
        thisWeight  = [];
        thisEstrous = [];
        thisChoice  = [];
        thisReward = [];
        thisSession = [];
        thisLatency_zscore = [];
        thisTrialTime = [];
        thisDates = [];
        thisEstrous_strict = [];
        sessionList = unique(session);
        thisTrial = []; 
        for nf = 1:numel(flist)
            try
                fid = fopen(fullfile(basefilename,ids{na}, flist{nf}));
                x=1;
                while x
                    l = fgetl(fid);
                    if strfind(l,'Start Date')
                        x=0;
                        thisDate = datenum(l(13:end));
                    end
                end
                fclose(fid);
            catch
                keyboard
            end
            if ismember(thisDate,dates)
                
                load(fullfile(basefilename,ids{na}, sprintf('%s_TAB.mat',flist{nf}(2:end-4))));
                sessionIdx  = find(session==sessionList(nf));
                if numel(data.trialIdx_all) ~= numel(sessionIdx)
                    keyboard
                end
                thisLatency = cat(1,thisLatency,trialStart(sessionIdx)');
                thisLatency_zscore = cat(1, thisLatency_zscore, nanzscore(trialStart(sessionIdx))');
                thisChoice  = cat(1,thisChoice,data.choice(data.trialIdx_all)');
                sessionIdx  = find(session==sessionList(nf));
                thisValue   = cat(1,thisValue,qChosenDiff(sessionIdx));
                thisValue_choice = cat(1,thisValue_choice,qDiff(sessionIdx)); 
                thisTrial   = cat(1,thisTrial,data.trialIdx_all'); 
                
                % consensus
                estrous = estrous_consensus(find(dates==thisDate),strcmp(estrous_consensus.Properties.VariableNames,ids{na}));
                estrous = eval(sprintf('estrous.%s;',ids{na}));
                if isempty(estrous)
                    estrous = nan(size(data.choice(data.trialIdx_all)'));
                else
                    estrous = (repmat(estrous,size(data.choice(data.trialIdx_all)')));
                end
                thisEstrous = cat(1,thisEstrous,estrous);
                
                % strict
                if isempty(find(dates_strict==thisDate))
                    estrous = (repmat({'NaN'},size(data.choice(data.trialIdx_all)')));
                else
                    
                    estrous = estrous_strict(find(dates_strict==thisDate),strcmp(estrous_strict.Properties.VariableNames,ids{na}));
                    
                    estrous = eval(sprintf('estrous.%s;',ids{na}));
                    if isempty(estrous)
                        estrous = nan(size(data.choice(data.trialIdx_all)'));
                    else
                        estrous = (repmat(estrous,size(data.choice(data.trialIdx_all)')));
                    end
                end
                thisEstrous_strict = cat(1,thisEstrous_strict,estrous);
                if numel(thisEstrous_strict) ~= numel(thisEstrous)
                    keyboard
                end
                % weight
                weight = weight_data(find(dates_weight==thisDate),strcmp(weight_data.Properties.VariableNames,ids{na}));
                weight = eval(sprintf('weight.%s;',ids{na}));
                thisWeight = cat(1,thisWeight,repmat(weight,size(data.choice(data.trialIdx_all)')));
                thisDates  = cat(1, thisDates, repmat(thisDate,size(data.choice(data.trialIdx_all)')));
                thisReward = cat(1,thisReward,data.reward(data.trialIdx_all)');
                thisSession = cat(1,thisSession,ones(size(data.reward(data.trialIdx_all)')).*nf);
                thisTrialTime = cat(1,thisTrialTime,data.trialStart(data.trialIdx_all)');
            end
        end
        
        Estrous = cat(1,Estrous,thisEstrous);
        Estrous_strict = cat(1,Estrous_strict,thisEstrous_strict);
        Weight = cat(1,Weight,thisWeight);
        Dates  = cat(1, Dates, thisDates);
        Value   = cat(1,Value,thisValue);
        Value_choice = cat(1,Value_choice,thisValue_choice);
        Choice  = cat(1,Choice, thisChoice);
        Latency = cat(1, Latency,thisLatency);
        Reward  = cat(1, Reward,thisReward);
        Latency_zscore = cat(1,Latency_zscore,thisLatency_zscore);
        ptiles = prctile(thisValue,linspace(0,100,binNum+1));
        TrialTime = cat(1,TrialTime, thisTrialTime);
        Trial = cat(1,Trial,thisTrial); 
        try
            [~,~,thisValue_ptile] = histcounts(thisValue,ptiles);
        catch
            keyboard
        end
        Value_ptile = cat(1, Value_ptile,thisValue_ptile);
        ptiles = prctile(thisValue_choice,linspace(0,100,binNum_choice+1));
        try
            [~,~,thisValue_ptile] = histcounts(thisValue_choice,ptiles);
        catch
            keyboard
        end
        Value_ptile_choice = cat(1, Value_ptile_choice,thisValue_ptile);
        Animal = cat(1,Animal,repmat(ids(na), size(thisValue)));
        Session = cat(1,Session,thisSession);
    end
    clear ids
    save(fullfile(whereAreWe('figurecode'),'processed_data','estrous',sprintf('estrousData_%s_%s',num2str(perfThresh),qFile)))
else
    load(fullfile(whereAreWe('figurecode'),'processed_data','estrous',sprintf('estrousData_%s_%s',num2str(perfThresh),qFile)))
end

Latency(Latency==0) = 0.01; % to make consistent with slower acquisition rate 


%% Number of sessions per animal 

for na = 1:numel(ids)
    thisDates = unique(Dates(strcmp(Animal,ids{na})));
    thisDates = thisDates(thisDates>=sessionStart);
    thisEstrous = arrayfun(@(x) Estrous(find(Dates==x&strcmp(Animal,ids{na}),1,'first')),thisDates);
    stageCount.E(na) = sum(strcmp(thisEstrous,'E'));
    stageCount.M(na) = sum(strcmp(thisEstrous,'M'));
    stageCount.D(na) = sum(strcmp(thisEstrous,'D'));
    stageCount.P(na) = sum(strcmp(thisEstrous,'P'));
end
stageCount.E_mean = mean(stageCount.E);
stageCount.E_std = std(stageCount.E);
stageCount.M_mean = mean(stageCount.M);
stageCount.M_std = std(stageCount.M);
stageCount.D_mean = mean(stageCount.D);
stageCount.D_std = std(stageCount.D);
stageCount.P_mean = mean(stageCount.P);
stageCount.P_std = std(stageCount.P);

%% Combine stages

Estrous_combo = cell((size(Estrous)));
Estrous_combo(contains(Estrous,'D')|contains(Estrous,'M')) = {'DM'};
Estrous_combo(contains(Estrous,'P')|contains(Estrous,'E')) = {'PE'};


%% Plot to threshold



Latency_thresh = Latency;

for na = 1:numel(ids)
    thisSession = unique(Session(strcmp(Animal,ids{na})));
    for ns = 1:numel(thisSession)
        sessionIdx = find(Session==thisSession(ns)&strcmp(Animal,ids{na})); 
        threshIdx = sessionIdx(find(Latency(sessionIdx)>intervalThresh,1,'first'));
        Latency_thresh(threshIdx:sessionIdx(end)) = NaN;
    end
end
%%
groups = {'PE';'DM'};
theseDates = sessionStart:stopDate;
clear latencyXestrous
for na = 1:numel(ids)
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
            for nb = 1:binNum
                eval(sprintf('latencyXestrous.%s(na,nb) = nanmean(Latency_thresh(ismember(Dates,theseDates)&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous,''%s'')|strcmp(Estrous,''%s''))));',gname,gname(1),gname(2)));
            end
        else
            gname = groups{ng};
            for nb = 1:binNum
                eval(sprintf('latencyXestrous.%s(na,nb) = nanmean(Latency_thresh(ismember(Dates,theseDates)&strcmp(Animal,ids{na})&Value_ptile==nb&strcmp(Estrous,''%s'')));',gname,groups{ng}));
            end
        end
    end
end


clear thisPlot
figure(); hold on
for ng = 1:numel(groups)
    if numel(groups{ng})>1
        gname = cat(2,groups{ng}(1), groups{ng}(end));
    else
        gname = groups{ng};
    end
    thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s',gname)));
    thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s(:,1)',gname)))));
end
b=bar(thisPlot');
pause(0.1)
for ng = 1:numel(groups)
    b(ng).EdgeColor = fmap(ng+1,:);
    b(ng).FaceColor = 'none';
    b(ng).LineWidth = 2;
    %  b(ng).FaceAlpha = .7;
end
for na = 1:numel(ids)
    clear aPlot xData bb
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
        else
            gname = groups{ng};
        end
        aPlot(ng,:) = eval(sprintf('latencyXestrous.%s(na,:)',gname));
        xData(ng,:) = b(ng).XData+b(ng).XOffset;
        bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng+1,:), 'MarkerEdgeColor','none');
    end
    
    p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
    uistack(bb,'top')
    
end

ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])
legend(groups)
%% Fit regression
X = table(categorical(Session), Latency_thresh, (Estrous), categorical(Value_ptile), (Animal),TrialTime,Dates,Weight,nanzscore(Trial),'VariableNames',{'Session','Latency','Estrous','Value','Animal','TrialTime','Dates','Weight','Trial'});
% X = X(strcmp(X.Estrous,'D')|strcmp(X.Estrous,'M')|strcmp(X.Estrous,'P')|strcmp(X.Estrous,'E')|strcmp(X.Estrous,'P/E')|strcmp(X.Estrous,'M/D')|strcmp(X.Estrous,'E/M')...
%     |strcmp(X.Estrous,'D/P'),:);
X = X(strcmp(X.Estrous,'D')|strcmp(X.Estrous,'M')|strcmp(X.Estrous,'P')|strcmp(X.Estrous,'E'),:);
X.Estrous(strcmp(X.Estrous,'D')|strcmp(X.Estrous,'M')) = {'DM'};
X.Estrous(strcmp(X.Estrous,'E')|strcmp(X.Estrous,'P')) = {'EP'};

X.Estrous = categorical(X.Estrous);
idx =cellfun(@(x) find(strcmp(X.Animal,x)),ids,'UniformOutput',false);
idx = cell2mat(idx);
X = X(idx,:);
X = X(X.Dates>=sessionStart,:);
X.Weight = nanzscore(X.Weight); 
X(isnan(X.Latency),:) = [];
X.Animal = categorical(X.Animal);


f = 'Latency ~ Value*Estrous  + (1+Value|Animal)';
stats.mdl = fitlme(X,f,'DummyVarCoding','effect');
stats.mdl_anova=dataset2struct(anova(stats.mdl,'DFMethod','Satterthwaite'));


% contrasts
h_DM = [1 1 0 0 0 0 0 0];
h_PE = [1 -1 0 0 0 0 0 0 ];


for nb = 1:binNum-1
    % DM v PE
    h1 = h_DM;
    idx = find(strcmp(stats.mdl.Coefficients.Name,sprintf('Value_%d',nb)));
    h1(idx) = 1;
    idx = find(contains(stats.mdl.Coefficients.Name,sprintf('Value_%d',nb))&contains(stats.mdl.Coefficients.Name,'DM'));
    h1(idx) = 1;
    h2 = h_PE;
    idx = find(strcmp(stats.mdl.Coefficients.Name,sprintf('Value_%d',nb)));
    h2(idx) = 1;
    idx = find(contains(stats.mdl.Coefficients.Name,sprintf('Value_%d',nb))&contains(stats.mdl.Coefficients.Name,'Estrous'));
    h2(idx) = -1;
    [stats.contrast.p(nb), stats.contrast.f(nb), stats.contrast.df1(nb),stats.contrast.df2(nb)] = coefTest(stats.mdl, h1-h2,0,'DFMethod','Satterthwaite');
end

% Value bin 4
nb = nb+1;
% DM v PE
h1 = h_DM;
idx = find(contains(stats.mdl.Coefficients.Name,'Value_'));
h1(idx(1:3)) = -1;
idx = find(contains(stats.mdl.Coefficients.Name,'Value_')&contains(stats.mdl.Coefficients.Name,'DM'));
h1(idx) = -1;
h2 = h_PE;
idx = find(contains(stats.mdl.Coefficients.Name,'Value_'));
h2(idx(1:3)) = -1;
idx = find(contains(stats.mdl.Coefficients.Name,'Value_')&contains(stats.mdl.Coefficients.Name,'Estrous_'));
h2(idx) = 1;
[stats.contrast.p(nb), stats.contrast.f(nb), stats.contrast.df1(nb),stats.contrast.df2(nb)] = coefTest(stats.mdl, h1-h2,0,'DFMethod','Satterthwaite');



end

