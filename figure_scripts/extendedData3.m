function stats=extendedData3

%% Parameters


perfThresh      = .1;
qFile           = 'qLearn_session_all.mat';
binNum          = 4;
ids = generateAnimalList('estrous');

intervalThresh =300;
sessionStart = datenum('3/1/22');
stopDate = datenum('5/15/2022');
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

[~,fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC);
fmap = fmap(round(linspace(20,size(fmap,1),3)),:);

%% Load data 

load(fullfile(whereAreWe('data'),sprintf('estrousData_%s_%s',num2str(perfThresh),qFile)))
Latency(Latency<0.01) = 0.01; % to make consistent with slower acquisition rate sessions

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
X = X(X.Dates<=stopDate,:);
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

