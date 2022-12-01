function stats=logLatencyHistograms_ctrl(behaviorTable,lowCut,highCut,ids_f,ids_m)
if lowCut < 0.01
    lowCut = 0.01;
end

behaviorTable=behaviorTable(~isnan(behaviorTable.trialInit_thresh),:);
%% Plot log normalized histograms by outcome and value
% Supp Fig 1C

% Plot log normalized histograms by outcome and value for each animal (control sessions)
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

% generate bins for histogram in log space
if highCut==inf
    bins = logspace(log10(lowCut),2,numel(log10(lowCut):2)*5);
else
    bins = logspace(log10(lowCut),log10(highCut),5*numel(-1:log10(highCut)));
end
bins =cat(2,bins,max(behaviorTable.trialInit));
% Previous reward female opsin

for na = 1:numel(ids_f)
    [femaleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==1&behaviorTable.previousReward==1&strcmp(behaviorTable.aID,ids_f{na}))),bins,'Normalization','probability');
end

figure();
subplot(2,2,1)
box off; hold on
set(gca,'XScale','log','YLim',[0 .2])
%plot mean
shadedErrorBar(bins(1:end-1), nanmean(femaleH,1), nanstd(femaleH)./sqrt(sum(~isnan(femaleH(:,1)))),'LineProps',{'Color',plotParams.femaleC, 'LineWidth',.25})
title('Reward')
%plot individuals
subplot(2,2,3)
box off; hold on
plot(repmat(bins(1:end-1)',1,size(femaleH,1)),femaleH','Color',[plotParams.femaleC .25], 'LineWidth',.25)
set(gca,'XScale','log','YLim',[0 .2])


% Previous no reward female opsin

for na = 1:numel(ids_f)
    [femaleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==1&behaviorTable.previousReward==0&strcmp(behaviorTable.aID,ids_f{na}))),bins,'Normalization','probability');
end
subplot(2,2,2)
box off; hold on
set(gca,'XScale','log','YLim',[0 .2])
shadedErrorBar(bins(1:end-1), nanmean(femaleH,1), nanstd(femaleH)./sqrt(sum(~isnan(femaleH(:,1)))),'LineProps',{'Color',plotParams.femaleC, 'LineWidth',.25})
title('No reward')
subplot(2,2,4)
box off; hold on
plot(repmat(bins(1:end-1)',1,size(femaleH,1)),femaleH','Color',[plotParams.femaleC .25], 'LineWidth',.25)
set(gca,'XScale','log','YLim',[0 .2])


% Previous reward male opsin

for na = 1:numel(ids_m)
    [maleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==0&strcmp(behaviorTable.aID,ids_m{na})&behaviorTable.previousReward==1)),bins,'Normalization','probability');
end
subplot(2,2,1)
box off; hold on
set(gca,'XScale','log','YLim',[0 .2])
shadedErrorBar(bins(1:end-1), nanmean(maleH,1), nanstd(maleH)./sqrt(sum(~isnan(maleH(:,1)))),'LineProps',{'Color',plotParams.maleC, 'LineWidth',1.5})
title('Reward')

% Previous no reward male opsin
for na = 1:numel(ids_m)
    [maleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==0&strcmp(behaviorTable.aID,ids_m{na})&behaviorTable.previousReward==0)),bins,'Normalization','probability');
end
subplot(2,2,3)
box off; hold on
plot(repmat(bins(1:end-1)',1,size(maleH,1)),maleH','Color',[plotParams.maleC .25], 'LineWidth',.5)
set(gca,'XScale','log','YLim',[0 .2])
subplot(2,2,2)
box off; hold on
set(gca,'XScale','log','YLim',[0 .2])
shadedErrorBar(bins(1:end-1), nanmean(maleH,1), nanstd(maleH)./sqrt(sum(~isnan(maleH(:,1)))),'LineProps',{'Color',plotParams.maleC, 'LineWidth',1.5})
title('No reward')
subplot(2,2,4)
box off; hold on
plot(repmat(bins(1:end-1)',1,size(maleH,1)),maleH','Color',[plotParams.maleC .25], 'LineWidth',.5)
set(gca,'XScale','log','YLim',[0 .2])

%% Plot in coarse bins
%Supp Fig S1D
clear femaleH maleH
bins= [0 1 10 max(behaviorTable.trialInit)];

vec = [];
groupvec = [];
% Previous reward

for na = 1:numel(ids_f)
    [femaleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==1&behaviorTable.previousReward==1&strcmp(behaviorTable.aID,ids_f{na}))),bins,'Normalization','probability');
end
for na = 1:numel(ids_m)
    [maleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==0&behaviorTable.previousReward==1&strcmp(behaviorTable.aID,ids_m{na}))),bins,'Normalization','probability');
end
figure();
subplot(2,2,1)
box off; hold on
%plot mean
mu = [nanmean(femaleH,1);nanmean(maleH,1)];
b=bar(mu');
b(1).EdgeColor = plotParams.femaleC;
b(2).EdgeColor = plotParams.maleC;
b(1).FaceColor = 'none';
b(2).FaceColor = 'none';
pause(0.01);
errorbar(b(1).XData+b(1).XOffset,nanmean(femaleH,1),nanstd(femaleH)./sqrt(sum(~isnan(femaleH(:,1)))),'Color',plotParams.femaleC,'LineStyle','none','CapSize',0)
errorbar(b(2).XData+b(2).XOffset,nanmean(maleH,1),nanstd(maleH)./sqrt(sum(~isnan(maleH(:,1)))),'Color',plotParams.maleC,'LineStyle','none','CapSize',0)
xaxis = repmat(b(1).XData+b(1).XOffset, size(femaleH,1),1);
plot(xaxis,femaleH,'o','MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerSize',4)
xaxis = repmat(b(2).XData+b(2).XOffset, size(maleH,1),1);
plot(xaxis,maleH,'o','MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerSize',4)
ylabel('Proportion of trials')
title('Reward')

for na = 1:size(femaleH,1)
    vec = cat(1,vec,femaleH(na,:)');
    groupvec = cat(1,groupvec,cat(2,ones(size(femaleH,2),1),[1:size(femaleH,2)]', ones(size(femaleH,2),1),ones(size(femaleH,2),1).*na));
end
acounter = na;
for na = 1:size(maleH,1)
    vec = cat(1,vec,maleH(na,:)');
    groupvec = cat(1,groupvec,cat(2,zeros(size(maleH,2),1),[1:size(maleH,2)]', ones(size(maleH,2),1),ones(size(maleH,2),1).*(na+acounter)));
end

for na = 1:numel(ids_f)
    [femaleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==1&behaviorTable.previousReward==0&strcmp(behaviorTable.aID,ids_f{na}))),bins,'Normalization','probability');
end
for na = 1:numel(ids_m)
    [maleH(na,:), ~] = histcounts((behaviorTable.trialInit_thresh(behaviorTable.female==0&behaviorTable.previousReward==0&strcmp(behaviorTable.aID,ids_m{na}))),bins,'Normalization','probability');
end
subplot(2,2,2)
box off; hold on
%plot mean
mu = [nanmean(femaleH,1);nanmean(maleH,1)];
b=bar(mu');
b(1).EdgeColor = plotParams.femaleC;
b(2).EdgeColor = plotParams.maleC;
b(1).FaceColor = 'none';
b(2).FaceColor = 'none';
pause(0.01);
errorbar(b(1).XData+b(1).XOffset,nanmean(femaleH,1),nanstd(femaleH)./sqrt(sum(~isnan(femaleH(:,1)))),'Color',plotParams.femaleC,'LineStyle','none','CapSize',0)
errorbar(b(2).XData+b(2).XOffset,nanmean(maleH,1),nanstd(maleH)./sqrt(sum(~isnan(maleH(:,1)))),'Color',plotParams.maleC,'LineStyle','none','CapSize',0)
xaxis = repmat(b(1).XData+b(1).XOffset, size(femaleH,1),1);
plot(xaxis,femaleH,'o','MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerSize',4)
xaxis = repmat(b(2).XData+b(2).XOffset, size(maleH,1),1);
plot(xaxis,maleH,'o','MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerSize',4)
ylabel('Proportion of trials')
title('No reward')

for na = 1:size(femaleH,1)
    vec = cat(1,vec,femaleH(na,:)');
    groupvec = cat(1,groupvec,cat(2,ones(size(femaleH,2),1),[1:size(femaleH,2)]', zeros(size(femaleH,2),1),ones(size(femaleH,2),1).*na));
end
acounter = na;
for na = 1:size(maleH,1)
    vec = cat(1,vec,maleH(na,:)');
    groupvec = cat(1,groupvec,cat(2,zeros(size(maleH,2),1),[1:size(maleH,2)]', zeros(size(maleH,2),1),ones(size(maleH,2),1).*(na+acounter)));
end


T = table(vec);
T.sex = (groupvec(:,1));
T.bin = (groupvec(:,2));
T.outcome = (groupvec(:,3));
T.animal = (groupvec(:,4));

[stats.p_short_rew,~,stats.stats_short_rew] = ranksum(T.vec(T.sex==1&T.bin==1&T.outcome==1),T.vec(T.sex==0&T.bin==1&T.outcome==1));
[stats.p_med_rew,~,stats.stats_med_rew] = ranksum(T.vec(T.sex==1&T.bin==2&T.outcome==1),T.vec(T.sex==0&T.bin==2&T.outcome==1));
[stats.p_long_rew,~,stats.stats_long_rew] = ranksum(T.vec(T.sex==1&T.bin==3&T.outcome==1),T.vec(T.sex==0&T.bin==3&T.outcome==1));
[stats.p_short_nrew,~,stats.stats_short_nrew] = ranksum(T.vec(T.sex==1&T.bin==1&T.outcome==0),T.vec(T.sex==0&T.bin==1&T.outcome==0));
[stats.p_med_nrew,~,stats.stats_med_nrew] = ranksum(T.vec(T.sex==1&T.bin==2&T.outcome==0),T.vec(T.sex==0&T.bin==2&T.outcome==0));
[stats.p_long_nrew,~,stats.stats_long_nrew] = ranksum(T.vec(T.sex==1&T.bin==3&T.outcome==0),T.vec(T.sex==0&T.bin==3&T.outcome==0));