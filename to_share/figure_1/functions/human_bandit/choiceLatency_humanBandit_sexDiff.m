function choiceLatency_humanBandit_sexDiff(cohort,zscoreFlag,cutoff,lowCutoff)

nback = 3; 

if nargin<3
    zscoreFlag = 0;
    cutoff = inf;
    lowCutoff = 0;
end


load(fullfile(whereAreWe('bucket'), 'Manuscript_figures','plotParams.mat'))

basefilename = fullfile(whereAreWe('bucket'),'human_bandit'); 

allHist_rew_f = [];
allHist_nrew_f = [];
allHist_rew_m = [];
allHist_nrew_m = []; 

% Convert data, if unconverted
dataConvert_humanBandit([cohort '_m']);
dataConvert_humanBandit([cohort '_f']);

histBins = linspace(0,2000,100); 

%% Female subjects
% Load each subject's data and calculate stay probability following rewarded and unrewarded outcomes 
sList = dir(fullfile(basefilename, [cohort '_f'], ['*.mat'])); 
sList = {sList(:).name};

for ns = 1:numel(sList)
   load(fullfile(basefilename, [cohort '_f'],sList{ns})); % load data trial
   saveFlag = 0; 
   if ~isfield(data, 'idx')
       data = getIndices_humanBandit(data); 
       saveFlag = 1;
   end
   
   if ~isfield(data,'multiRew')
       data = multiOutcome_humanBandit(data);
       saveFlag = 1;
   end
   
   if zscoreFlag == 3
       [rewCDF_f(ns,:),histBins] = histcounts(log(data.choiceLatency(data.idx.prevRew)), 100,'Normalization','cdf');
       allHist_rew_f = cat(1,allHist_rew_f,log(data.choiceLatency(data.idx.prevRew)));
       nrewCDF_f(ns,:) = histcounts(log(data.choiceLatency(data.idx.prevNRew)), histBins,'Normalization','cdf');
       allHist_nrew_f = cat(1,allHist_nrew_f,log(data.choiceLatency(data.idx.prevNRew)));
   else
       rewCDF_f(ns,:) = histcounts(data.choiceLatency(data.idx.prevRew), histBins,'Normalization','cdf');
       allHist_rew_f = cat(1,allHist_rew_f,data.choiceLatency(data.idx.prevRew));
       nrewCDF_f(ns,:) = histcounts(data.choiceLatency(data.idx.prevNRew), histBins,'Normalization','cdf');
       allHist_nrew_f = cat(1,allHist_nrew_f,data.choiceLatency(data.idx.prevNRew));
   end
   
   data.choiceLatency(data.choiceLatency>cutoff |data.choiceLatency<lowCutoff) = NaN; 
   data.choiceLatency_zscore(data.choiceLatency>cutoff |data.choiceLatency<lowCutoff) = NaN; 
   if zscoreFlag == 1
       data.choiceLatency_zscore = nanzscore(data.choiceLatency);
       choiceLatency_prevRew(ns) = nanmean(data.choiceLatency_zscore(data.idx.prevRew));
       choiceLatency_prevNRew(ns) = nanmean(data.choiceLatency_zscore(data.idx.prevNRew));
       for nb = 1:nback
           choiceLatency_multiRew(ns,nb) = nanmean(data.choiceLatency_zscore(data.multiRew{nb}));
           choiceLatency_multiNRew(ns,nb) = nanmean(data.choiceLatency_zscore(data.multiNRew{nb}));
       end
       std_rew(ns) = nanstd(data.choiceLatency_zscore(data.idx.prevRew))./sqrt(numel(data.idx.prevRew));
       std_nrew(ns) = nanstd(data.choiceLatency_zscore(data.idx.prevNRew))./sqrt(numel(data.idx.prevNRew));
   elseif zscoreFlag == 0
       choiceLatency_prevRew(ns) = nanmean(data.choiceLatency(data.idx.prevRew));
       choiceLatency_prevNRew(ns) = nanmean(data.choiceLatency(data.idx.prevNRew));
       for nb = 1:nback
           choiceLatency_multiRew(ns,nb) = nanmean(data.choiceLatency(data.multiRew{nb}));
           choiceLatency_multiNRew(ns,nb) = nanmean(data.choiceLatency(data.multiNRew{nb}));
       end
       std_rew(ns) = nanstd(data.choiceLatency(data.idx.prevRew))./sqrt(numel(data.idx.prevRew));
       std_nrew(ns) = nanstd(data.choiceLatency(data.idx.prevNRew))./sqrt(numel(data.idx.prevNRew));
   elseif zscoreFlag == 2
       choiceLatency_prevRew(ns) = nanmean(data.choiceLatency(data.idx.prevRew)-nanmean(data.choiceLatency));
       choiceLatency_prevNRew(ns) = nanmean(data.choiceLatency(data.idx.prevNRew)-nanmean(data.choiceLatency));
       for nb = 1:nback
           choiceLatency_multiRew(ns,nb) = nanmean(data.choiceLatency(data.multiRew{nb})-nanmean(data.choiceLatency));
           choiceLatency_multiNRew(ns,nb) = nanmean(data.choiceLatency(data.multiNRew{nb})-nanmean(data.choiceLatency));
       end
       std_rew(ns) = nanstd(data.choiceLatency(data.idx.prevRew)-nanmean(data.choiceLatency))./sqrt(numel(data.idx.prevRew));
       std_nrew(ns) = nanstd(data.choiceLatency(data.idx.prevNRew)-nanmean(data.choiceLatency))./sqrt(numel(data.idx.prevNRew));
   elseif zscoreFlag == 3
       choiceLatency_prevRew(ns) = nanmean(log(data.choiceLatency(data.idx.prevRew)));
       choiceLatency_prevNRew(ns) = nanmean(log(data.choiceLatency(data.idx.prevNRew)));
       for nb = 1:nback
           choiceLatency_multiRew(ns,nb) = nanmean(log(data.choiceLatency(data.multiRew{nb})));
           choiceLatency_multiNRew(ns,nb) = nanmean(log(data.choiceLatency(data.multiNRew{nb})));
       end
       std_rew(ns) = nanstd(log(data.choiceLatency(data.idx.prevRew)))./sqrt(numel(data.idx.prevRew));
       std_nrew(ns) = nanstd(log(data.choiceLatency(data.idx.prevNRew)))./sqrt(numel(data.idx.prevNRew));
   end
   
   
%    
%    if saveFlag
%        save(fullfile(basefilename,sList{ns}),'data');
%    end
end

%% male subjects
% Load each subject's data and calculate stay probability following rewarded and unrewarded outcomes 
sList = dir(fullfile(basefilename, [cohort '_m'], ['*.mat'])); 
sList = {sList(:).name};


for ns = 1:numel(sList)
   load(fullfile(basefilename, [cohort '_m'],sList{ns})); % load data trial
   saveFlag = 0; 
   if ~isfield(data, 'idx')
       data = getIndices_humanBandit(data); 
       saveFlag = 1;
   end
   
   if ~isfield(data,'multiRew')
       data = multiOutcome_humanBandit(data);
       saveFlag = 1;
   end

    if zscoreFlag == 3
       [rewCDF_m(ns,:),histBins] = histcounts(log(data.choiceLatency(data.idx.prevRew)), 100,'Normalization','cdf');
       allHist_rew_m = cat(1,allHist_rew_m,log(data.choiceLatency(data.idx.prevRew)));
       nrewCDF_m(ns,:) = histcounts(log(data.choiceLatency(data.idx.prevNRew)), histBins,'Normalization','cdf');
       allHist_nrew_m = cat(1,allHist_nrew_f,log(data.choiceLatency(data.idx.prevNRew)));
   else
       rewCDF_m(ns,:) = histcounts(data.choiceLatency(data.idx.prevRew), histBins,'Normalization','cdf');
       allHist_rew_m = cat(1,allHist_rew_m,data.choiceLatency(data.idx.prevRew));
       nrewCDF_m(ns,:) = histcounts(data.choiceLatency(data.idx.prevNRew), histBins,'Normalization','cdf');
       allHist_nrew_m = cat(1,allHist_nrew_m,data.choiceLatency(data.idx.prevNRew));
   end
   
   data.choiceLatency(data.choiceLatency>cutoff | data.choiceLatency<lowCutoff) = NaN;
   
   if zscoreFlag == 1
       data.choiceLatency_zscore = nanzscore(data.choiceLatency);
       choiceLatency_prevRew_m(ns) = nanmean(data.choiceLatency_zscore(data.idx.prevRew));
       choiceLatency_prevNRew_m(ns) = nanmean(data.choiceLatency_zscore(data.idx.prevNRew));
       for nb = 1:nback
           choiceLatency_multiRew_m(ns,nb) = nanmean(data.choiceLatency_zscore(data.multiRew{nb}));
           choiceLatency_multiNRew_m(ns,nb) = nanmean(data.choiceLatency_zscore(data.multiNRew{nb}));
       end
       
       std_rew_m(ns) = nanstd(data.choiceLatency_zscore(data.idx.prevRew))./sqrt(numel(data.idx.prevRew));
       std_nrew_m(ns) = nanstd(data.choiceLatency_zscore(data.idx.prevNRew))./sqrt(numel(data.idx.prevNRew));
       
   elseif zscoreFlag == 0
       choiceLatency_prevRew_m(ns) = nanmean(data.choiceLatency(data.idx.prevRew));
       choiceLatency_prevNRew_m(ns) = nanmean(data.choiceLatency(data.idx.prevNRew));
       try
       for nb = 1:nback
           choiceLatency_multiRew_m(ns,nb) = nanmean(data.choiceLatency(data.multiRew{nb}));
           try
           choiceLatency_multiNRew_m(ns,nb) = nanmean(data.choiceLatency(data.multiNRew{nb}));
           catch
           end
       end
       catch
           keyboard
       end
       std_rew_m(ns) = nanstd(data.choiceLatency(data.idx.prevRew))./sqrt(numel(data.idx.prevRew));
       std_nrew_m(ns) = nanstd(data.choiceLatency(data.idx.prevNRew))./sqrt(numel(data.idx.prevNRew));
   elseif zscoreFlag == 2
       choiceLatency_prevRew_m(ns) = nanmean(data.choiceLatency(data.idx.prevRew)-nanmean(data.choiceLatency));
       choiceLatency_prevNRew_m(ns) = nanmean(data.choiceLatency(data.idx.prevNRew)-nanmean(data.choiceLatency));
       for nb = 1:nback
           choiceLatency_multiRew_m(ns,nb) = nanmean(data.choiceLatency(data.multiRew{nb})-nanmean(data.choiceLatency));
           choiceLatency_multiNRew_m(ns,nb) = nanmean(data.choiceLatency(data.multiNRew{nb})-nanmean(data.choiceLatency));
       end
       std_rew_m(ns) = nanstd(data.choiceLatency(data.idx.prevRew)-nanmean(data.choiceLatency))./sqrt(numel(data.idx.prevRew));
       std_nrew_m(ns) = nanstd(data.choiceLatency(data.idx.prevNRew)-nanmean(data.choiceLatency))./sqrt(numel(data.idx.prevNRew));
   elseif zscoreFlag == 3
       choiceLatency_prevRew_m(ns) = nanmean(log(data.choiceLatency(data.idx.prevRew)));
       choiceLatency_prevNRew_m(ns) = nanmean(log(data.choiceLatency(data.idx.prevNRew)));
       for nb = 1:nback
           choiceLatency_multiRew_m(ns,nb) = nanmean(log(data.choiceLatency(data.multiRew{nb})));
           choiceLatency_multiNRew_m(ns,nb) = nanmean(log(data.choiceLatency(data.multiNRew{nb})));
       end
       std_rew_m(ns) = nanstd(log(data.choiceLatency(data.idx.prevRew)))./sqrt(numel(data.idx.prevRew));
       std_nrew_m(ns) = nanstd(log(data.choiceLatency(data.idx.prevNRew)))./sqrt(numel(data.idx.prevNRew));
   end
   
   
   
%    
%    if saveFlag
%        save(fullfile(basefilename,sList{ns}),'data');
%    end
end

%% Plot distribution of male and female data 
if zscoreFlag == 3
    [rew_f] = histcounts(allHist_rew_f,histBins,'Normalization','cdf');
    nrew_f = histcounts(allHist_nrew_f,histBins,'Normalization','cdf');
    rew_m = histcounts(allHist_rew_m,histBins,'Normalization','cdf');
    nrew_m = histcounts(allHist_nrew_m,histBins,'Normalization','cdf');
else
    rew_f = histcounts(allHist_rew_f,histBins,'Normalization','cdf');
    nrew_f = histcounts(allHist_nrew_f,histBins,'Normalization','cdf');
    rew_m = histcounts(allHist_rew_m,histBins,'Normalization','cdf');
    nrew_m = histcounts(allHist_nrew_m,histBins,'Normalization','cdf');
end

figure(); hold on

plot(histBins(1:end-1),rew_f,'Color',femaleC);
plot(histBins(1:end-1),nrew_f,':','Color',femaleC);
plot(histBins(1:end-1),rew_m,'Color',maleC);
plot(histBins(1:end-1),nrew_m,':','Color',maleC);
xlabel('Choice latency (ms)')
ylabel('Proportion of trials')
legend({'Female: prev rew';'Female: prev unrew'; 'Male: prev rew';'Male: prev unrew'},'Box','off'); 




figure(); hold on
p(1)=shadedErrorBar(histBins(1:end-1), mean(rewCDF_f), std(rewCDF_f)./sqrt(size(rewCDF_f,1)),'LineProps',{'Color',femaleC});
p(2)=shadedErrorBar(histBins(1:end-1), mean(nrewCDF_f), std(nrewCDF_f)./sqrt(size(nrewCDF_f,1)),'LineProps',{':','Color',femaleC});

p(3)=shadedErrorBar(histBins(1:end-1), mean(rewCDF_m), std(rewCDF_m)./sqrt(size(rewCDF_m,1)),'LineProps',{'Color',maleC});
p(4)=shadedErrorBar(histBins(1:end-1), mean(nrewCDF_m), std(nrewCDF_m)./sqrt(size(nrewCDF_m,1)),'LineProps',{':','Color',maleC});

xlabel('Choice latency (ms)')
ylabel('Proportion of trials')
legend([p(1).mainLine p(2).mainLine p(3).mainLine p(4).mainLine],{'Female: prev rew';'Female: prev unrew'; 'Male: prev rew';'Male: prev unrew'},'Box','off'); 


%% Plot choice latency by previous outcome

% figure();
% g=[ones(size(choiceLatency_prevRew)) ones(size(choiceLatency_prevRew_m)).*2  ones(size(choiceLatency_prevNRew)).*3 ones(size(choiceLatency_prevNRew_m)).*4];
% gc=[ones(size(choiceLatency_prevRew)) ones(size(choiceLatency_prevRew_m)).*2  ones(size(choiceLatency_prevNRew)).*1 ones(size(choiceLatency_prevNRew_m)).*2];
% 
% b=boxplot([choiceLatency_prevRew choiceLatency_prevRew_m choiceLatency_prevNRew choiceLatency_prevNRew_m],g,'PlotStyle','compact',...
%     'ColorGroup',gc,'Jitter',0);
% box off
% legend({'Male';'Female'})
% % 
% % figure()
% % b = boxchart(g,[choiceLatency_prevRew choiceLatency_prevRew_m choiceLatency_prevNRew choiceLatency_prevNRew_m])
% 
% data.prevRew_f = choiceLatency_prevRew;
% data.prevNRew_f = choiceLatency_prevNRew;
% data.prevRew_m = choiceLatency_prevRew_m;
% data.prevNRew_m = choiceLatency_prevNRew_m; 
% figure(), v=violinplot(data,{'Rew female';'Rew male';'NRew female';'NRew male'},'GroupOrder',{'prevRew_f';'prevRew_m';'prevNRew_f';'prevNRew_m'},'Bandwidth',50,'ViolinColor',femaleC);
% 
% 
 mu = [nanmean(choiceLatency_prevRew) nanmean(choiceLatency_prevNRew); nanmean(choiceLatency_prevRew_m) nanmean(choiceLatency_prevNRew_m)]; 
 sem = [std(choiceLatency_prevRew )./sqrt(numel(choiceLatency_prevRew)) std(choiceLatency_prevNRew)./sqrt(numel(choiceLatency_prevNRew)); std(choiceLatency_prevRew_m)./sqrt(numel(choiceLatency_prevRew_m)) std(choiceLatency_prevNRew_m)./sqrt(numel(choiceLatency_prevNRew_m))];

figure()
b=bar(mu');
b(1).FaceColor = 'none';
b(1).EdgeColor = femaleC;
b(1).LineWidth = 2;
b(2).FaceColor = 'none';
b(2).EdgeColor = maleC;
b(2).LineWidth = 2;
box off
hold on
pause(.01)
errorbar(b(1).XData+b(1).XOffset, b(1).YData, sem(1,:),'LineStyle','none','Color',femaleC,'CapSize',0,'LineWidth',2)
errorbar(b(2).XData+b(2).XOffset, b(2).YData, sem(2,:),'LineStyle','none','Color',maleC,'CapSize',0,'LineWidth',2)

plot(repmat((b(1).XData(1)+b(1).XOffset)',1,numel(choiceLatency_prevRew))+randi([-4 4],1,numel(choiceLatency_prevRew)).*.01, [choiceLatency_prevRew],'o','MarkerFaceColor',femaleC,'MarkerEdgeColor',femaleC,'MarkerSize',5)
plot(repmat((b(1).XData(2)+b(1).XOffset)',1,numel(choiceLatency_prevRew))+randi([-4 4],1,numel(choiceLatency_prevRew)).*.01, [choiceLatency_prevNRew],'o','MarkerFaceColor',femaleC,'MarkerEdgeColor',femaleC,'MarkerSize',5)

plot(repmat((b(2).XData(1)+b(2).XOffset)',1,numel(choiceLatency_prevRew_m))+randi([-4 4],1,numel(choiceLatency_prevRew_m)).*.01, [choiceLatency_prevRew_m],'o','MarkerFaceColor',maleC,'MarkerEdgeColor',maleC,'MarkerSize',5)
plot(repmat((b(2).XData(2)+b(2).XOffset)',1,numel(choiceLatency_prevRew_m))+randi([-4 4],1,numel(choiceLatency_prevRew_m)).*.01, [choiceLatency_prevNRew_m],'o','MarkerFaceColor',maleC,'MarkerEdgeColor',maleC,'MarkerSize',5)


set(gca, 'FontSize', 14, 'XTickLabel', {'reward'; 'no reward'})
if zscoreFlag
    ylabel('Choice latency (zscore)','FontSize',16)
else
    ylabel('Choice latency (ms)','FontSize',16)
end
xlabel('Previous trial outcome','FontSize',16)
legend({'female';'male'})

% 2-way ANOVA
% concatenate data
vec = cat(1,choiceLatency_prevRew',choiceLatency_prevNRew',choiceLatency_prevRew_m',choiceLatency_prevNRew_m'); %concatenate female previously rewarded, unrewarded, male previously rewarded unrewarded latencies by subject
% grouping variables
groupvec = cat(1,cat(2,ones(size(choiceLatency_prevRew')),ones(size(choiceLatency_prevRew'))),... % female = 1, previous reward = 1
    cat(2,ones(size(choiceLatency_prevNRew')),zeros(size(choiceLatency_prevNRew'))),... % female = 1, previous reward = 0
    cat(2,zeros(size(choiceLatency_prevRew_m')),ones(size(choiceLatency_prevRew_m'))),... % female = 0, previous reward = 1
    cat(2,zeros(size(choiceLatency_prevNRew_m')),zeros(size(choiceLatency_prevNRew_m')))); % female = 0, previous reward = 1
[stats.p_anova,stats.t_anova,stats.s_anova]= anovan(vec, groupvec,'varnames', {'sex';'outcome'}, 'display', 'off', 'model','full');
% post hoc test 
figure();stats.postHoc_tukey = multcompare(stats.s_anova,'dimension',[1 2],'CType','hsd','Display','on');
%%
figure();
if ~zscoreFlag
    hold on, plot([200 700],[200, 700],':k')
    title('Choice latency (ms)')
else
    hold on
    plot([0 0], [-.3 .3], ':k')
    plot([-.3 .3], [0 0], ':k')
%plot([-.4 .4],[-.4,.4],':k')
    title('Choice latency (zscore)')
end
hold on
scatter(choiceLatency_prevRew,choiceLatency_prevNRew,30,femaleC,'filled')
plot([choiceLatency_prevRew; choiceLatency_prevRew], [std_nrew+choiceLatency_prevNRew; -std_nrew+choiceLatency_prevNRew], 'Color',femaleC)
plot([std_rew+choiceLatency_prevRew; -std_rew+choiceLatency_prevRew],[choiceLatency_prevNRew; choiceLatency_prevNRew], 'Color',femaleC)

scatter(choiceLatency_prevRew_m,choiceLatency_prevNRew_m,30,maleC,'filled')
plot([choiceLatency_prevRew_m; choiceLatency_prevRew_m], [std_nrew_m+choiceLatency_prevNRew_m; -std_nrew_m+choiceLatency_prevNRew_m], 'Color',maleC)
plot([std_rew_m+choiceLatency_prevRew_m; -std_rew_m+choiceLatency_prevRew_m],[choiceLatency_prevNRew_m; choiceLatency_prevNRew_m],  'Color',maleC)
set(gca, 'FontSize', 14)
xlabel('Previous reward','FontSize',16)
ylabel('Previous no reward','FontSize',16)
box off

%% Plot Choice latency by multiple outcomes 

figure()
hold on

errorbar(1:nback, nanmean(choiceLatency_multiRew,1), std(choiceLatency_multiRew,[],1)./sqrt(numel(sList)),'Color', femaleC, 'CapSize',0)
errorbar(1:nback, nanmean(choiceLatency_multiNRew,1),std(choiceLatency_multiNRew,[],1)./sqrt(numel(sList)),'Color',femaleC,'LineStyle', ':', 'CapSize',0)

errorbar(1:nback, nanmean(choiceLatency_multiRew_m,1), std(choiceLatency_multiRew_m,[],1)./sqrt(numel(sList)),'Color', maleC, 'CapSize',0)
errorbar(1:nback, nanmean(choiceLatency_multiNRew_m,1),std(choiceLatency_multiNRew_m,[],1)./sqrt(numel(sList)),'Color',maleC,'LineStyle', ':', 'CapSize',0)

set(gca, 'XTick', 1:nback,'XLim', [.5 nback+.5],'FontSize',14)
ylabel('Choice latency (ms)', 'FontSize',16)
xlabel('Number of previous outcomes','FontSize',16)
legend({'female: rewarded';'female: unrewarded';'male: rewarded';'male: unrewarded'}); 
%% Plot difference between rewarded and unrewarded latencies

mu = [nanmean(choiceLatency_prevNRew-choiceLatency_prevRew); nanmean(choiceLatency_prevNRew_m-choiceLatency_prevRew_m)]; 
 sem = [std(choiceLatency_prevNRew-choiceLatency_prevRew )./sqrt(numel(choiceLatency_prevRew)) ; std(choiceLatency_prevNRew_m-choiceLatency_prevRew_m)./sqrt(numel(choiceLatency_prevRew_m))];

figure(); hold on
b=bar(mu(1));
b(2)= bar(2,mu(2));
b(1).FaceColor = 'none';
b(1).EdgeColor = femaleC;
b(1).LineWidth = 2;
b(2).FaceColor = 'none';
b(2).EdgeColor = maleC;
b(2).LineWidth = 2;
box off
hold on
pause(.01)
errorbar(b(1).XData+b(1).XOffset, b(1).YData, sem(1,:),'LineStyle','none','Color',femaleC,'CapSize',0,'LineWidth',2)
errorbar(b(2).XData+b(2).XOffset, b(2).YData, sem(2,:),'LineStyle','none','Color',maleC,'CapSize',0,'LineWidth',2)

plot(repmat((b(1).XData(1)+b(1).XOffset)',1,numel(choiceLatency_prevRew))+randi([-4 4],1,numel(choiceLatency_prevRew)).*.01, [choiceLatency_prevNRew-choiceLatency_prevRew],'o','MarkerFaceColor',femaleC,'MarkerEdgeColor',femaleC,'MarkerSize',5)

plot(repmat((b(2).XData(1)+b(2).XOffset)',1,numel(choiceLatency_prevRew_m))+randi([-4 4],1,numel(choiceLatency_prevRew_m)).*.01, [choiceLatency_prevNRew_m-choiceLatency_prevRew_m],'o','MarkerFaceColor',maleC,'MarkerEdgeColor',maleC,'MarkerSize',5)


set(gca, 'FontSize', 14)%, 'XTickLabel', {'reward'; 'no reward'})
if zscoreFlag
    ylabel('Previous no reward - previous reward (zscore)','FontSize',16)
else
    ylabel('Previous no reward - previous reward (ms)','FontSize',16)
end
set(gca,'XTick',[1 2],'XTickLabel',{'Female';'Male'})
