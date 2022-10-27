%% Estrous histogram

triggerState = 'D';

estBins = [-1:5];
latencyMean = nan(numel(ids),numel(estBins));
propD = nan(numel(ids),numel(estBins));
propE = nan(numel(ids),numel(estBins));
propM = nan(numel(ids),numel(estBins));
propP = nan(numel(ids),numel(estBins));
sessionStop = datenum('5/15/2022');
sessionStart = datenum('3/1/2022');
latencyTrials = [];
valueTrials = [];
intervalTrials = [];
animalTrials = [];
estrousTrials = [];
for na = 1:numel(ids)
    clear thisLatencyMean thisInterval thisDay
    sessionEstrous    =  eval(sprintf('estrous_consensus.%s;',ids{na}));
    thisDates   = unique(datenum(estrous_consensus.Date));
    sessionEstrous(thisDates<sessionStart) = [];
    thisDates(thisDates<sessionStart) = [];
    sessionEstrous(thisDates>sessionStop) = [];
    thisDates(thisDates>sessionStop) = [];
    
    idx = find(strcmp(sessionEstrous,triggerState));
    idx(find(diff(idx)==1)) = [];
    thisDates = thisDates(idx); 
    thisInterval = nan(size(thisDates));
    thisLatencyMean = nan(size(thisDates,1),numel(estBins)); 
    thisDay = nan(size(thisDates,1),numel(estBins)); 
    for ns = 1:numel(idx)
        x=1;
        for nb = estBins
            thisIdx = find(strcmp(Animal,ids{na})&Dates==thisDates(ns)+nb);
            if ~isempty(thisIdx)
            thisLatencyMean(ns,x) = nanmean(Latency_thresh(thisIdx));
            if strcmp(unique(Estrous(thisIdx)),'P')
                thisDay(ns,x) = 1;
            elseif strcmp(unique(Estrous(thisIdx)),'E')
                thisDay(ns,x) = 2;
            elseif strcmp(unique(Estrous(thisIdx)),'M')
                thisDay(ns,x) = 3;
             elseif strcmp(unique(Estrous(thisIdx)),'D')
                thisDay(ns,x) = 4;
            else
                thisDay(ns,x) = NaN;      
            end
            latencyTrials = cat(1,latencyTrials,Latency_thresh(thisIdx));
            valueTrials = cat(1,valueTrials,Value(thisIdx));
            animalTrials = cat(1,animalTrials,ones(size(thisIdx)).*na);
            intervalTrials = cat(1,intervalTrials,ones(size(thisIdx)).*nb); 
            estrousTrials  = cat(1,estrousTrials,ones(size(thisIdx)).*thisDay(ns,x));
            end
            
            x=x+1;
            
        end
    end
    latencyMean(na,:) = nanmean(thisLatencyMean);
    propP(na,:) = sum(thisDay==1,1)./sum(~isnan(thisDay),1);
    propE(na,:) = sum(thisDay==2,1)./sum(~isnan(thisDay),1);
    propM(na,:) = sum(thisDay==3,1)./sum(~isnan(thisDay),1);
    propD(na,:) = sum(thisDay==4,1)./sum(~isnan(thisDay),1);
end
propE(propE==Inf) = NaN;
propM(propM==Inf) = NaN;
propP(propP==Inf) = NaN;
propD(propD==Inf) = NaN;



f='latency~value+(1+value|animal)';

clear value value_se intercept intercept_se
x=1;
for nb = (estBins)
   % for ng = 1;
   for ng = 1:max(estrousTrials)
        thisIdx = find(intervalTrials==nb&estrousTrials==ng);
       % thisIdx = find(intervalTrials==nb);
        X = table(latencyTrials(thisIdx),valueTrials(thisIdx),animalTrials(thisIdx), 'VariableNames',{'latency','value','animal'});
        try
            mdl = fitlme(X,f);
            value(ng,x) = mdl.Coefficients.Estimate(2);
            value_se(ng,x) = mdl.Coefficients.SE(2);
            intercept(ng,x) = mdl.Coefficients.Estimate(1);
            intercept_se(ng,x) = mdl.Coefficients.SE(1);
        catch
            value(ng,x) = NaN;
            value_se(ng,x) = NaN;
            intercept(ng,x) = NaN;
            intercept_se(ng,x) = NaN;
        end
    end
    x=x+1;
    
end


figure(), subplot(3,2,1); hold on 
errorbar(estBins(1:end),nanmean(latencyMean),nanstd(latencyMean)./sqrt(sum(~isnan(latencyMean),1)),'k','LineWidth',1,'CapSize',0,'LineStyle','none')
%plot(estBins(1:end-1),nanmean(latencyMean),'o','MarkerFaceColor','k','MarkerEdgeColor','none')
X = table(reshape(latencyMean,numel(latencyMean),1),'VariableNames',{'latency'});
X.interval = categorical(reshape(repmat(estBins(1:end),numel(ids),1),numel(latencyMean),1));
X.Animal = repmat([1:numel(ids)]',numel(estBins),1);
mdl = fitlme(X,'latency~interval+(1|Animal)');
anova(mdl)
% mdl = fitlm(X,'latency~interval');
% p=plot(mdl);
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = 'k';
% p(1).MarkerEdgeColor = 'none';
% p(1).MarkerSize = 3;
% p(2).Color = 'k';
% p(3).Color = 'k';
% p(4).Color = 'k';

for na = 1:size(latencyMean,1)
scatter(estBins(1:end),latencyMean(na,:), 15,'MarkerFaceColor', 'k','MarkerEdgeColor','none','MarkerFaceAlpha',.5)
end
    
xlabel('Days since estrous')
ylabel('Average trial initiation latency (s)')
set(gca,'XLim',[estBins(1)-1 estBins(end)+1]);
legend('off')

subplot(3,2,5)
hold on
errorbar(estBins(1:end), nanmean(propP), nanstd(propP)./sqrt(sum(~isnan(propP))), 'Color', fmap(1,:),'CapSize',0,'LineWidth',1)
errorbar(estBins(1:end), nanmean(propE), nanstd(propE)./sqrt(sum(~isnan(propE))), 'Color', fmap(2,:),'CapSize',0,'LineWidth',1)
errorbar(estBins(1:end), nanmean(propM), nanstd(propM)./sqrt(sum(~isnan(propM))), 'Color', fmap(3,:),'CapSize',0,'LineWidth',1)

% errorbar(estBins(1:end-1), nanmean(propM+propD), nanstd(propM+propD)./sqrt(sum(~isnan(propM+propD))), 'Color', fmap(2,:),'CapSize',0,'LineWidth',1)
% errorbar(estBins(1:end-1), nanmean(propE+propP), nanstd(propE+propP)./sqrt(sum(~isnan(propE+propP))), 'Color', fmap(4,:),'CapSize',0,'LineWidth',1)

errorbar(estBins(1:end), nanmean(propD), nanstd(propD)./sqrt(sum(~isnan(propD))), 'Color', fmap(4,:),'CapSize',0,'LineWidth',1)
legend({'P';'E';'M';'D'})
%legend({'D/M'  ;'P/E'})
ylabel('Proportion of sessions')
xlabel('Days since estrous')
set(gca,'XLim',[estBins(1)-1 estBins(end)+1]);
box off
size(value)


subplot(3,2,3); hold on
for ng = 1:max(estrousTrials)
    errorbar(estBins(1:end),value(ng,:),value_se(ng,:),'Color',fmap(ng,:),'CapSize',0,'LineWidth',1)
end
%errorbar(estBins(1:end),value,value_se,'Color','k','CapSize',0,'LineWidth',1)
ylabel('Value coefficient')
xlabel('Days since estrous')
set(gca,'XLim',[estBins(1)-1 estBins(end)+1]);

% errorbar(estBins(1:end),value,value_se,'Color','k','CapSize',0,'LineWidth',1,'LineStyle','none')
% plot(estBins(1:end),value,'ok','MarkerFaceColor','k','MarkerEdgeColor','none')
% ylabel('Value coefficient')
% xlabel('Days since estrous')
% set(gca,'XLim',[estBins(1)-1 estBins(end)+1]);
% box off
% X = table(value','VariableNames',{'value'});
% X.interval = categorical(estBins(1:end)');
% mdl = fitlm(X,'value~interval');
% anova(mdl)
% p=plot(mdl);
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = 'k';
% p(1).MarkerEdgeColor = 'none';
% p(1).MarkerSize = 3;
% p(2).Color = 'k';
% p(3).Color = 'k';
% p(4).Color = 'k';
% legend('off')

% figure();
% subplot(1,2,1)
% errorbar(estBins(1:end-1),value,value_se,'Color','k','CapSize',0,'LineWidth',1)
% ylabel('Value coefficient')
% xlabel('Days since estrous')
% set(gca,'XLim',[-1 estBins(end-1)+1]);
% subplot(1,2,2)
% errorbar(estBins(1:end-1),intercept,intercept_se,'Color','k','CapSize',0,'LineWidth',1)
% ylabel('Intercept')
% xlabel('Days since estrous')
% set(gca,'XLim',[-1 estBins(end-1)+1]);


%% Estrous interval

estBins = [0:15,50];
latencyMean = nan(numel(ids),numel(estBins)-1);
propD = nan(numel(ids),numel(estBins)-1);
propE = nan(numel(ids),numel(estBins)-1);
propM = nan(numel(ids),numel(estBins)-1);
propP = nan(numel(ids),numel(estBins)-1);
sessionStop = datenum('5/15/2022');
sessionStart = datenum('3/1/2022');
latencyTrials = [];
valueTrials = [];
intervalTrials = [];
animalTrials = [];
for na = 1:numel(ids)
    clear thisLatencyMean thisInterval thisDay
    sessionEstrous    =  eval(sprintf('estrous_consensus.%s;',ids{na}));
    thisDates   = unique(datenum(estrous_consensus.Date));
    idx            = find(cellfun(@(x) isempty(x),  eval(sprintf('estrous_consensus.%s;',ids{na}))));     
    thisDates(idx) = [];
    sessionEstrous(idx) = [];
    sessionEstrous(thisDates<sessionStart) = [];
    thisDates(thisDates<sessionStart) = [];
        sessionEstrous(thisDates>sessionStop) = [];
    thisDates(thisDates>sessionStop) = [];
    
    breaks = diff(thisDates);
    breaks = cat(1,1,find(breaks>1)+1,numel(thisDates));
    thisInterval{na} = nan(size(thisDates));
    thisLatencyMean{na} = nan(size(thisDates)); 
    thisDay = cell(size(thisDates));
    for nb = 2:numel(breaks)
        firstEstrous   = find(strcmp(sessionEstrous(breaks(nb-1):breaks(nb)),'E'),1,'first')+breaks(nb-1)-1;
        for ns = firstEstrous:breaks(nb)-1
            thisIdx = find(strcmp(Animal,ids{na})&Dates==thisDates(ns));
            if ~isempty(thisIdx)
                
            prevEstrous = find(strcmp(sessionEstrous(1:ns),'E'),1,'last');
            thisLatencyMean{na}(ns) = nanmean(Latency_thresh(thisIdx));
            thisInterval{na}(ns) = ns-prevEstrous;
            latencyTrials = cat(1,latencyTrials,Latency_thresh(thisIdx));
            valueTrials = cat(1,valueTrials,Value(thisIdx));
            thisDay(ns) = unique(Estrous(thisIdx));
            intervalTrials = cat(1,intervalTrials,ones(size(thisIdx)).*thisInterval{na}(ns));
            animalTrials = cat(1,animalTrials,ones(size(thisIdx)).*na);
            end
        end
    end
   
    [~,~,thisBins] = histcounts(thisInterval{na},estBins);
    
    for nb = 1:numel(estBins)-1
       latencyMean(na,nb) = nanmean(thisLatencyMean{na}(thisBins==nb));
       binIdx = find(thisBins==nb);
       binIdx = binIdx(ismember(binIdx,idx));
       propD(na,nb) = sum(strcmp(thisDay(thisBins==nb),'D'))/sum(thisBins==nb); 
       propE(na,nb) = sum(strcmp(thisDay(thisBins==nb),'E'))/sum(thisBins==nb); 
       propM(na,nb) = sum(strcmp(thisDay(thisBins==nb),'M'))/sum(thisBins==nb);  
       propP(na,nb) = sum(strcmp(thisDay(thisBins==nb),'P'))/sum(thisBins==nb); 
    end
end
[~,~,intervalTrials2] = histcounts(intervalTrials,estBins); 


f='latency~value+(1+value|animal)';


clear value value_se intercept intercept_se
for nb = 1:numel(estBins)-1
    thisIdx = find(intervalTrials2==nb);
    X = table(latencyTrials(thisIdx),valueTrials(thisIdx),animalTrials(thisIdx), 'VariableNames',{'latency','value','animal'});
    try
        mdl = fitlme(X,f);
        value(nb) = mdl.Coefficients.Estimate(2);
        value_se(nb) = mdl.Coefficients.SE(2);
        intercept(nb) = mdl.Coefficients.Estimate(1);
        intercept_se(nb) = mdl.Coefficients.SE(1);
    catch
        value(nb) = NaN;
        value_se(nb) = NaN;
        intercept(nb) = NaN;
        intercept_se(nb) = NaN;
    end
end


figure(), subplot(3,2,1); hold on 
errorbar(estBins(1:end-1),nanmean(latencyMean),nanstd(latencyMean)./sqrt(sum(~isnan(latencyMean),1)),'k','LineWidth',1,'CapSize',0,'LineStyle','none')
%plot(estBins(1:end-1),nanmean(latencyMean),'o','MarkerFaceColor','k','MarkerEdgeColor','none')
X = table(reshape(latencyMean,numel(latencyMean),1),'VariableNames',{'latency'});
X.interval = reshape(repmat(estBins(1:end-1),numel(ids),1),numel(latencyMean),1);
X.Animal = repmat([1:numel(ids)]',numel(estBins)-1,1);
mdl = fitlme(X,'latency~interval+(1|Animal)') 
X = table(reshape(latencyMean,numel(latencyMean),1),'VariableNames',{'latency'});
X.interval = reshape(repmat(estBins(1:end-1),numel(ids),1),numel(latencyMean),1);
X.Animal = repmat([1:numel(ids)]',numel(estBins)-1,1);
mdl = fitlm(X,'latency~interval');
p=plot(mdl);
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'k';
p(1).MarkerEdgeColor = 'none';
p(1).MarkerSize = 3;
p(2).Color = 'k';
p(3).Color = 'k';
p(4).Color = 'k';

% for na = 1:size(latencyMean,1)
% scatter(estBins(1:end-1),latencyMean(na,:), 15,'MarkerFaceColor', 'k','MarkerEdgeColor','none','MarkerFaceAlpha',.5)
% end
    
xlabel('Days since estrous')
ylabel('Average trial initiation latency (s)')
set(gca,'XLim',[-1 estBins(end-1)+1]);
legend('off')

subplot(3,2,5)
hold on
errorbar(estBins(1:end-1), nanmean(propP), nanstd(propP)./sqrt(sum(~isnan(propP))), 'Color', fmap(1,:),'CapSize',0,'LineWidth',1)
errorbar(estBins(1:end-1), nanmean(propE), nanstd(propE)./sqrt(sum(~isnan(propE))), 'Color', fmap(2,:),'CapSize',0,'LineWidth',1)
errorbar(estBins(1:end-1), nanmean(propM), nanstd(propM)./sqrt(sum(~isnan(propM))), 'Color', fmap(3,:),'CapSize',0,'LineWidth',1)

% errorbar(estBins(1:end-1), nanmean(propM+propD), nanstd(propM+propD)./sqrt(sum(~isnan(propM+propD))), 'Color', fmap(2,:),'CapSize',0,'LineWidth',1)
% errorbar(estBins(1:end-1), nanmean(propE+propP), nanstd(propE+propP)./sqrt(sum(~isnan(propE+propP))), 'Color', fmap(4,:),'CapSize',0,'LineWidth',1)

errorbar(estBins(1:end-1), nanmean(propD), nanstd(propD)./sqrt(sum(~isnan(propD))), 'Color', fmap(4,:),'CapSize',0,'LineWidth',1)
legend({'P';'E';'M';'D'})
%legend({'D/M';'P/E'})
ylabel('Proportion of sessions')
xlabel('Days since estrous')
set(gca,'XLim',[-1 estBins(end-1)+1]);
box off



subplot(3,2,3); hold on
errorbar(estBins(1:end-1),value,value_se,'Color','k','CapSize',0,'LineWidth',1,'LineStyle','none')
plot(estBins(1:end-1),value,'ok','MarkerFaceColor','k','MarkerEdgeColor','none')
ylabel('Value coefficient')
xlabel('Days since estrous')
set(gca,'XLim',[-1 estBins(end-1)+1]);
box off
X = table(value','VariableNames',{'value'});
X.interval = estBins(1:end-1)';
mdl = fitlm(X,'value~interval');
p=plot(mdl);
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'k';
p(1).MarkerEdgeColor = 'none';
p(1).MarkerSize = 3;
p(2).Color = 'k';
p(3).Color = 'k';
p(4).Color = 'k';
legend('off')
% figure();
% subplot(1,2,1)
% errorbar(estBins(1:end-1),value,value_se,'Color','k','CapSize',0,'LineWidth',1)
% ylabel('Value coefficient')
% xlabel('Days since estrous')
% set(gca,'XLim',[-1 estBins(end-1)+1]);
% subplot(1,2,2)
% errorbar(estBins(1:end-1),intercept,intercept_se,'Color','k','CapSize',0,'LineWidth',1)
% ylabel('Intercept')
% xlabel('Days since estrous')
% set(gca,'XLim',[-1 estBins(end-1)+1]);

%% Inter estrous interval


estBins = [1:2:16,50];
latencyMean = nan(numel(ids),numel(estBins)-1);
latencyTrials = [];
valueTrials = [];
intervalTrials = [];
animalTrials = [];
for na = 1:numel(ids)
    clear thisLatencyMean thisInterval
    thisIdx        = find(strcmp(Animal,ids{na})&Dates>=sessionStart);
    sessionEstrous    =  eval(sprintf('estrous_consensus.%s;',ids{na}));
    thisDates   = unique(datenum(estrous_consensus.Date));
    idx            = find(cellfun(@(x) isempty(x),  eval(sprintf('estrous_consensus.%s;',ids{na}))));     
    thisDates(idx) = [];
    sessionEstrous(idx) = [];
    breaks = diff(thisDates);
    breaks = cat(1,1,find(breaks>1)+1,thisDates(end));
    thisInterval = nan(size(thisDates));
    thisLatencyMean = nan(size(thisDates)); 
    
    for nb = 2:numel(breaks)
                firstEstrous   = find(strcmp(sessionEstrous(breaks(nb-1):breaks(nb)),'E'),1,'first')+breaks(nb-1)-1;
                lastEstrous   = find(strcmp(sessionEstrous(breaks(nb-1):breaks(nb)),'E'),1,'last')+breaks(nb-1)-1;

    for ns = firstEstrous:lastEstrous
        thisIdx = find(strcmp(Animal,ids{na})&Dates==thisDates(ns));
        prevEstrous = find(strcmp(sessionEstrous(1:ns),'E'),1,'last');
        nextEstrous = find(strcmp(sessionEstrous(ns:end),'E'),1,'first');
        thisLatencyMean(ns) = nanmean(Latency_thresh(thisIdx));
        thisInterval(ns) = nextEstrous-1+ns-prevEstrous;
        thisDay(ns) = thisDates(ns); 
        latencyTrials = cat(1,latencyTrials,Latency_thresh(thisIdx));
        valueTrials = cat(1,valueTrials,Value(thisIdx));
        intervalTrials = cat(1,intervalTrials,ones(size(thisIdx)).*thisInterval(ns));
        animalTrials = cat(1,animalTrials,ones(size(thisIdx)).*na);
    end
    end
   
    thisLatencyMean(1:firstEstrous-1) = [];
    thisInterval(1:firstEstrous-1) = [];
    
    [~,~,thisBins] = histcounts(thisInterval,estBins);
    for nb = 1:numel(estBins)-1
       latencyMean(na,nb) = nanmean(thisLatencyMean(thisBins==nb)); 
    end

end
figure(), errorbar(estBins(1:end-1),nanmean(latencyMean),nanstd(latencyMean)./sqrt(sum(~isnan(latencyMean),1)),'k','LineWidth',1)
xlabel('inter-estrous interval (days)')
ylabel('Average trial initiation latency (s)')
set(gca,'XLim',[-1 estBins(end-1)+1]);
f='latency~value+(1+value|animal)';
clear value value_se intercept intercept_se
for nb = 1:numel(estBins)-1
    thisIdx = find(intervalTrials==nb);
    try
        X = table(latencyTrials(thisIdx),valueTrials(thisIdx),animalTrials(thisIdx), 'VariableNames',{'latency','value','animal'});
        mdl = fitlme(X,f);
        value(nb) = mdl.Coefficients.Estimate(2);
        value_se(nb) = mdl.Coefficients.SE(2);
        intercept(nb) = mdl.Coefficients.Estimate(1);
        intercept_se(nb) = mdl.Coefficients.SE(1);
    catch
        value(nb) = nan;
        value_se(nb) = nan;
        intercept(nb) = nan;
        intercept_se(nb) = nan;
    end
end

figure();
subplot(1,2,1)
errorbar(estBins(1:end-1),value,value_se,'Color','k','CapSize',0,'LineWidth',1)
ylabel('Value coefficient')
xlabel('inter-estrous interval (days)')
set(gca,'XLim',[-1 estBins(end-1)+1]);
subplot(1,2,2)
errorbar(estBins(1:end-1),intercept,intercept_se,'Color','k','CapSize',0,'LineWidth',1)
ylabel('Intercept')
xlabel('inter-estrous interval (days)')
set(gca,'XLim',[-1 estBins(end-1)+1]);

%% state by state regression
groups = {'D';'P';'E';'M'};
%groups = {'PE';'DM'};
f = 'Latency ~ Value + (1+Value|Animal) + (1+Value|Session)';
sessionStart = datenum('3/1/2022');
sessionStop = datenum('5/15/2022');
for ng = 1:numel(groups)
    if numel(groups{ng})>1
        idx = find(Dates>=sessionStart&Dates<=sessionStop&(strcmp(Estrous,groups{ng}(1))|strcmp(Estrous,groups{ng}(2))));
    else
        idx = find(Dates>=sessionStart&Dates<=sessionStop&strcmp(Estrous,groups{ng}));
    end
    X = table(Latency_thresh(idx), Value(idx), (Animal(idx)),categorical(Session(idx)), 'VariableNames',{'Latency';'Value';'Animal';'Session'});
    idx =cellfun(@(x) find(strcmp(X.Animal,x)),ids,'UniformOutput',false);
    idx = cell2mat(idx);
    X = X(idx,:);
    mdl = fitlme(X,f);
    fits.value(ng) = mdl.Coefficients.Estimate(2);
    fits.value_se(ng) = mdl.Coefficients.SE(2);
    fits.intercept(ng) = mdl.Coefficients.Estimate(1);
    fits.intercept_se(ng) = mdl.Coefficients.SE(1);
    
    [re, renames] = randomEffects(mdl);
    thisIds = unique(renames.Level(strcmp(renames.Group,'Animal')));
    for na = 1:numel(thisIds)
        fits.re_value{ng}(na) = re(strcmp(renames.Level,thisIds{na})&strcmp(renames.Name,'Value')&strcmp(renames.Group,'Animal'))+fits.value(ng);
        fits.re_intercept{ng}(na) = re(strcmp(renames.Level,thisIds{na})&strcmp(renames.Name,'(Intercept)')&strcmp(renames.Group,'Animal'))+fits.intercept(ng);
    end
end

figure(); 
subplot(1,2,1); hold on
for ng = 1:numel(groups)
   errorbar(ng,fits.value(ng),fits.value_se(ng),'o','Color',fmap(ng,:),'MarkerFaceColor',fmap(ng,:),'MarkerEdgeColor','none','LineWidth',1)
    
end
set(gca,'XLim',[0 numel(groups)+1],'XTick',1:numel(groups),'XTickLabel',groups)
xlabel('Estrous stage')
ylabel('Value coefficient')
subplot(1,2,2); hold on
for ng = 1:numel(groups)
   errorbar(ng,fits.intercept(ng),fits.intercept_se(ng),'o','Color',fmap(ng,:),'MarkerFaceColor',fmap(ng,:),'MarkerEdgeColor','none','LineWidth',1)
    
end
set(gca,'XLim',[0 numel(groups)+1],'XTick',1:numel(groups),'XTickLabel',groups)
xlabel('Estrous stage')
ylabel('Intercept')

figure();
subplot(1,2,1); hold on
for na = 1:numel(ids)
for ng = 1:numel(groups)
   scatter(ones(size(fits.re_value{ng})).*ng, fits.re_value{ng},20,'MarkerFaceColor',fmap(ng,:),'MarkerEdgeColor','none')    
   aplot(ng) = fits.re_value{ng}(na);
end
b=plot(1:numel(groups), aplot,'Color',[.7 .7 .7]);
uistack(b,'bottom')
end
set(gca,'XLim',[0 numel(groups)+1],'XTick',1:numel(groups),'XTickLabel',groups)
xlabel('Estrous stage')
ylabel('Value coefficient')

subplot(1,2,2); hold on
for na = 1:numel(ids)
for ng = 1:numel(groups)
   scatter(ones(size(fits.re_intercept{ng})).*ng, fits.re_intercept{ng},20,'MarkerFaceColor',fmap(ng,:),'MarkerEdgeColor','none')    
   aplot(ng) = fits.re_intercept{ng}(na);
end
b=plot(1:numel(groups), aplot,'Color',[.7 .7 .7]);
uistack(b,'bottom')
end
set(gca,'XLim',[0 numel(groups)+1],'XTick',1:numel(groups),'XTickLabel',groups)
xlabel('Estrous stage')
ylabel('Intercept')


% %% Plot by time in session not by value 
% groups = {'D';'P';'E';'M'};
% clear latencyXestrous 
% 
% timeBins = linspace(0,5400,2);
% [~,~,TimeBins] = histcounts(TrialTime,timeBins);
% 
% for na = 1:numel(ids)
%     for ntb = 1:numel(timeBins)-1
%         
%         
%             for ng = 1:numel(groups)
%                 if numel(groups{ng})>1
%                     gname = cat(2,groups{ng}(1), groups{ng}(end));
%                 else
%                     gname = groups{ng};
%                 end
%                 eval(sprintf('latencyXestrous.%s_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&strcmp(Estrous,''%s'')));',gname,groups{ng}));
%             end
%             latencyXestrous.est_pro_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&(strcmp(Estrous,'P')|strcmp(Estrous,'E'))));
%             latencyXestrous.di_met_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&(strcmp(Estrous,'M')|strcmp(Estrous,'D'))));
%             
%             latencyXestrous.overall_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})));
%             
%             
%         end
%     
% end
% clear thisPlot
% figure();
% for nt = 1:numel(timeBins)-1
%         subplot(1,numel(timeBins)-1,nt); hold on
% 
%     for ng = 1:numel(groups)
%         if numel(groups{ng})>1
%             gname = cat(2,groups{ng}(1), groups{ng}(end));
%         else
%             gname = groups{ng};
%         end
%         thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s_time{nt}',gname)));
%         thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s_time{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s_time{nt}(:,1)',gname)))));
%         b(ng) = bar(ng,thisPlot(ng,:));
%     end
%    % b=bar(thisPlot');
%     pause(0.1)
%     for ng = 1:numel(groups)
%         b(ng).EdgeColor = fmap(ng,:);
%         b(ng).FaceColor = 'none';
%         b(ng).LineWidth = 2;
%       %  b(ng).FaceAlpha = .7;
%     end
%     for na = 1:numel(ids)
%         clear aPlot xData
%         for ng = 1:numel(groups)
%             if numel(groups{ng})>1
%                 gname = cat(2,groups{ng}(1), groups{ng}(end));
%             else
%                 gname = groups{ng};
%             end
%             aPlot(ng,:) = eval(sprintf('latencyXestrous.%s_time{nt}(na)',gname));
%             xData(ng,:) = b(ng).XData+b(ng).XOffset;
%             bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
%         end
%         
%        p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
%        uistack(bb,'top') 
%         
%     end
% end
% ylabel('Trial initiation latency (s)')
% set(gca,'XTick',[1:binNum])
% legend(groups)
% 
% for ng = 1:numel(groups)
%     for ngg = 1:numel(groups)
%         if ng~=ngg
%            eval(sprintf('[stats.p_%s_%s, ~, stats.stats_%s_%s] = signrank(latencyXestrous.%s_time{1},latencyXestrous.%s_time{1});',groups{ng},groups{ngg}, groups{ng},groups{ngg},groups{ng},groups{ngg}));
%         end
%     end
% end


%% only include certain sessions


%% Plot by time in session not by value 
% groups = {'D';'P';'E';'M'};
% clear latencyXestrous 
% theseDates = datenum(['4/25/2022';'4/26/2022';'5/10/2022';'5/11/2022';'5/12/2022';'5/13/2022';'5/14/2022']);
% 
% timeBins = linspace(0,5400,2);
% [~,~,TimeBins] = histcounts(TrialTime,timeBins);
% 
% for na = 1:numel(ids)
%     for ntb = 1:numel(timeBins)-1
%         
%         
%             for ng = 1:numel(groups)
%                 if numel(groups{ng})>1
%                     gname = cat(2,groups{ng}(1), groups{ng}(end));
%                 else
%                     gname = groups{ng};
%                 end
%                 eval(sprintf('latencyXestrous.%s_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&ismember(Dates,theseDates)&strcmp(Animal,ids{na})&strcmp(Estrous,''%s'')));',gname,groups{ng}));
%             end
%             latencyXestrous.est_pro_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&ismember(Dates,theseDates)&strcmp(Animal,ids{na})&(strcmp(Estrous,'P')|strcmp(Estrous,'E'))));
%             latencyXestrous.di_met_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&ismember(Dates,theseDates)&strcmp(Animal,ids{na})&(strcmp(Estrous,'M')|strcmp(Estrous,'D'))));
%             
%             latencyXestrous.overall_time{ntb}(na) = nanmean(Latency(TimeBins==ntb&ismember(Dates,theseDates)&strcmp(Animal,ids{na})));
%             
%             
%         end
%     
% end
% clear thisPlot
% figure();
% for nt = 1:numel(timeBins)-1
%         subplot(1,numel(timeBins)-1,nt); hold on
% 
%     for ng = 1:numel(groups)
%         if numel(groups{ng})>1
%             gname = cat(2,groups{ng}(1), groups{ng}(end));
%         else
%             gname = groups{ng};
%         end
%         thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s_time{nt}',gname)));
%         thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s_time{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s_time{nt}(:,1)',gname)))));
%         b(ng) = bar(ng,thisPlot(ng,:));
%     end
%    % b=bar(thisPlot');
%     pause(0.1)
%     for ng = 1:numel(groups)
%         b(ng).EdgeColor = fmap(ng,:);
%         b(ng).FaceColor = 'none';
%         b(ng).LineWidth = 2;
%       %  b(ng).FaceAlpha = .7;
%     end
%     for na = 1:numel(ids)
%         clear aPlot xData
%         for ng = 1:numel(groups)
%             if numel(groups{ng})>1
%                 gname = cat(2,groups{ng}(1), groups{ng}(end));
%             else
%                 gname = groups{ng};
%             end
%             aPlot(ng,:) = eval(sprintf('latencyXestrous.%s_time{nt}(na)',gname));
%             xData(ng,:) = b(ng).XData+b(ng).XOffset;
%             bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
%         end
%         
%        p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
%        uistack(bb,'top') 
%         
%     end
% end
% ylabel('Trial initiation latency (s)')
% set(gca,'XTick',[1:binNum])
% legend(groups)
% 
% for ng = 1:numel(groups)
%     for ngg = 1:numel(groups)
%         if ng~=ngg
%            eval(sprintf('[stats.p_%s_%s, ~, stats.stats_%s_%s] = signrank(latencyXestrous.%s_time{1},latencyXestrous.%s_time{1});',groups{ng},groups{ngg}, groups{ng},groups{ngg},groups{ng},groups{ngg}));
%         end
%     end
% end
% 


%% Plot w/o value zscored 
% clear latencyXestrous 
% Latency_Z = nan(size(Latency));
% timeThresh = 3600;
% for na = 1:numel(ids)
%       idx = find(strcmp(Animal,ids{na})&TrialTime<=timeThresh&Dates>=sessionStart);
%       Latency_Z(idx) = nanzscore(Latency(idx));
% end
% 
% groups = {'D';'P';'E';'M'};
% clear latencyXestrous 
% 
% timeBins = linspace(0,3600,2);
% [~,~,TimeBins] = histcounts(TrialTime,timeBins);
% 
% for na = 1:numel(ids)
%     for ntb = 1:numel(timeBins)-1
%         
%         
%             for ng = 1:numel(groups)
%                 if numel(groups{ng})>1
%                     gname = cat(2,groups{ng}(1), groups{ng}(end));
%                 else
%                     gname = groups{ng};
%                 end
%                 eval(sprintf('latencyXestrous.%s_time{ntb}(na) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&strcmp(Estrous,''%s'')));',gname,groups{ng}));
%             end
%             latencyXestrous.est_pro_time{ntb}(na) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&(strcmp(Estrous,'P')|strcmp(Estrous,'E'))));
%             latencyXestrous.di_met_time{ntb}(na) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&(strcmp(Estrous,'M')|strcmp(Estrous,'D'))));
%             
%             latencyXestrous.overall_time{ntb}(na) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})));
%             
%             
%         end
%     
% end
% clear thisPlot
% figure();
% for nt = 1:numel(timeBins)-1
%         subplot(1,numel(timeBins)-1,nt); hold on
% 
%     for ng = 1:numel(groups)
%         if numel(groups{ng})>1
%             gname = cat(2,groups{ng}(1), groups{ng}(end));
%         else
%             gname = groups{ng};
%         end
%         thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s_time{nt}',gname)));
%         thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s_time{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s_time{nt}(:,1)',gname)))));
%         b(ng) = bar(ng,thisPlot(ng,:));
%     end
%    % b=bar(thisPlot');
%     pause(0.1)
%     for ng = 1:numel(groups)
%         b(ng).EdgeColor = fmap(ng,:);
%         b(ng).FaceColor = 'none';
%         b(ng).LineWidth = 2;
%       %  b(ng).FaceAlpha = .7;
%     end
%     for na = 1:numel(ids)
%         clear aPlot xData
%         for ng = 1:numel(groups)
%             if numel(groups{ng})>1
%                 gname = cat(2,groups{ng}(1), groups{ng}(end));
%             else
%                 gname = groups{ng};
%             end
%             aPlot(ng,:) = eval(sprintf('latencyXestrous.%s_time{nt}(na)',gname));
%             xData(ng,:) = b(ng).XData+b(ng).XOffset;
%             bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
%         end
%         
%        p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
%        uistack(bb,'top') 
%         
%     end
% end
% ylabel('Trial initiation latency (z-score)')
% set(gca,'XTick',[1:binNum])
% legend(groups)

%% Plot by time in session 
groups = {'D';'P';'E';'M'};
clear latencyXestrous 

timeBins = linspace(0,7200,3);
[~,~,TimeBins] = histcounts(TrialTime,timeBins);

for na = 1:numel(ids)
    for ntb = 1:numel(timeBins)-1
        
        for nb = 1:binNum
            for ng = 1:numel(groups)
                if numel(groups{ng})>1
                    gname = cat(2,groups{ng}(1), groups{ng}(end));
                else
                    gname = groups{ng};
                end
                eval(sprintf('latencyXestrous.%s_time{ntb}(na,nb) = nanmean(Latency_thresh(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&strcmp(Estrous,''%s'')));',gname,groups{ng}));
            end
            latencyXestrous.est_pro_time{ntb}(na,nb) = nanmean(Latency_thresh(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous,'P')|strcmp(Estrous,'E'))));
            latencyXestrous.di_met_time{ntb}(na,nb) = nanmean(Latency_thresh(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous,'M')|strcmp(Estrous,'D'))));
            
            latencyXestrous.overall_time{ntb}(na,nb) = nanmean(Latency_thresh(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb));
            
            
        end
    end
end
clear thisPlot thisError
figure();
for nt = 1:numel(timeBins)-1
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
        else
            gname = groups{ng};
        end
        thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s_time{nt}',gname)));
        thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s_time{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s_time{nt}(:,1)',gname)))));
    end
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(0.1)
    for ng = 1:numel(groups)
        b(ng).EdgeColor = fmap(ng,:);
        b(ng).FaceColor = 'none';
        b(ng).LineWidth = 2;
      %  b(ng).FaceAlpha = .7;
    end
    for na = 1:numel(ids)
        clear aPlot xData
        for ng = 1:numel(groups)
            if numel(groups{ng})>1
                gname = cat(2,groups{ng}(1), groups{ng}(end));
            else
                gname = groups{ng};
            end
            aPlot(ng,:) = eval(sprintf('latencyXestrous.%s_time{nt}(na,:)',gname));
            xData(ng,:) = b(ng).XData+b(ng).XOffset;
            bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
        end
        
       p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
       uistack(bb,'top') 
        
    end
end
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])
legend(groups)

clear thisPlot aPlot xData
figure();
for nt = 1:numel(timeBins)-1
    
    thisPlot(2,:) = nanmean(latencyXestrous.est_pro_time{nt});
    thisPlot(1,:) = nanmean(latencyXestrous.di_met_time{nt});
   
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(.1)
    for ng = 1:2
        b(ng).EdgeColor = b(ng).FaceColor;
        b(ng).FaceColor = 'none';
    end
    for na = 1:numel(ids)
        aPlot(2,:) = latencyXestrous.est_pro_time{nt}(na,:);
        aPlot(1,:) = latencyXestrous.di_met_time{nt}(na,:);
        xData(2,:) = b(2).XData+b(2).XOffset;
        xData(1,:) = b(1).XData+b(1).XOffset;


        
       plot(xData,aPlot,'o-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.25)
        clear aPlot xData
        
    end
end
legend(b,{'di/met';'est/pro'})
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])



%% Plot by time in session strict
groups = {'D';'P';'E';'M'};
%groups = {'D';'D/P';'P';'P/E';'E';'E/M';'M';'M/D'};
[~,fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC);
fmap = fmap(round(linspace(20,size(fmap,1),numel(groups))),:);

clear latencyXestrous 

timeBins = linspace(0,7200,2);
[~,~,TimeBins] = histcounts(TrialTime,timeBins);

for na = 1:numel(ids)
    for ntb = 1:numel(timeBins)-1
        
        for nb = 1:binNum
            for ng = 1:numel(groups)
                if numel(groups{ng})>1
                    gname = cat(2,groups{ng}(1), groups{ng}(end));
                else
                    gname = groups{ng};
                end
                eval(sprintf('latencyXestrous.%s_time{ntb}(na,nb) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&strcmp(Estrous_strict,''%s'')));',gname,groups{ng}));
            end
            latencyXestrous.est_pro_time{ntb}(na,nb) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous_strict,'P')|strcmp(Estrous_strict,'E'))));
            latencyXestrous.di_met_time{ntb}(na,nb) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous_strict,'M')|strcmp(Estrous_strict,'D'))));
            
            latencyXestrous.overall_time{ntb}(na,nb) = nanmean(Latency(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb));
            
            
        end
    end
end
clear thisPlot
figure();
for nt = 1:numel(timeBins)-1
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
        else
            gname = groups{ng};
        end
        thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s_time{nt}',gname)));
        thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s_time{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s_time{nt}(:,1)',gname)))));
    end
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(0.1)
    for ng = 1:numel(groups)
        b(ng).EdgeColor = fmap(ng,:);
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
            aPlot(ng,:) = eval(sprintf('latencyXestrous.%s_time{nt}(na,:)',gname));
            xData(ng,:) = b(ng).XData+b(ng).XOffset;
            bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
        end
        
       p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
       uistack(bb,'top') 
        
    end
end
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])
legend(groups)

clear thisPlot
figure();
for nt = 1:numel(timeBins)-1
    
    thisPlot(2,:) = nanmean(latencyXestrous.est_pro_time{nt});
    thisPlot(1,:) = nanmean(latencyXestrous.di_met_time{nt});
   
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(.1)
    for ng = 1:2
        b(ng).EdgeColor = b(ng).FaceColor;
        b(ng).FaceColor = 'none';
    end
    for na = 1:numel(ids)
        aPlot(2,:) = latencyXestrous.est_pro_time{nt}(na,:);
        aPlot(1,:) = latencyXestrous.di_met_time{nt}(na,:);
        xData(2,:) = b(2).XData+b(2).XOffset;
        xData(1,:) = b(1).XData+b(1).XOffset;


        
       plot(xData,aPlot,'o-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.25)
        clear aPlot xData
        
    end
end
legend(b,{'di/met';'est/pro'})
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])








%% Plot by time in session, z-score by animal 
clear latencyXestrous 
Latency_Z = nan(size(Latency));
timeThresh = 3600;
for na = 1:numel(ids)
      idx = find(strcmp(Animal,ids{na})&TrialTime<=timeThresh&Dates>=sessionStart);
      Latency_Z(idx) = nanzscore(Latency(idx));
end

groups = {'D';'P';'E';'M'};

[~,fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC);
fmap = fmap(round(linspace(20,size(fmap,1),numel(groups))),:);

timeBins = linspace(0,3600,2);
[~,~,TimeBins] = histcounts(TrialTime,timeBins);
for na = 1:numel(ids)
    for ntb = 1:numel(timeBins)-1
        
        for nb = 1:binNum
            for ng = 1:numel(groups)
                if numel(groups{ng})>1
                    gname = cat(2,groups{ng}(1), groups{ng}(end));
                else
                    gname = groups{ng};
                end
                eval(sprintf('latencyXestrous.%s{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&strcmp(Estrous,''%s'')));',gname,groups{ng}));
            end
            latencyXestrous.est_pro_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous,'P')|strcmp(Estrous,'E'))));
            latencyXestrous.di_met_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous,'M')|strcmp(Estrous,'D'))));
            
            latencyXestrous.overall_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb));
            
        end
    end
end
clear thisPlot
figure();
clear bb
for nt = 1:numel(timeBins)-1
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
        else
            gname = groups{ng};
        end
        thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s{nt}',gname)));
        thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s{nt}(:,1)',gname)))));
    end
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(0.1)
    for ng = 1:numel(groups)
        b(ng).EdgeColor = fmap(ng,:);
        b(ng).FaceColor = 'none';
        b(ng).LineWidth = 2;
      %  b(ng).FaceAlpha = .7;
    end
    for na = 1:numel(ids)
        clear aPlot xData
        for ng = 1:numel(groups)
            if numel(groups{ng})>1
                gname = cat(2,groups{ng}(1), groups{ng}(end));
            else
                gname = groups{ng};
            end
            aPlot(ng,:) = eval(sprintf('latencyXestrous.%s{nt}(na,:)',gname));
            xData(ng,:) = b(ng).XData+b(ng).XOffset;
            bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
        end
        
       p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
       uistack(bb,'top') 
        
    end
end
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])
legend(groups)

clear thisPlot
figure();
for nt = 1:numel(timeBins)-1
    
    thisPlot(2,:) = nanmean(latencyXestrous.est_pro_time{nt});
    thisPlot(1,:) = nanmean(latencyXestrous.di_met_time{nt});
   
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    for ng = 1:2
        b(ng).EdgeColor = b(ng).FaceColor;
        b(ng).FaceColor = 'none';
    end
    for na = 1:numel(ids)
        aPlot(2,:) = latencyXestrous.est_pro_time{nt}(na,:);
        aPlot(1,:) = latencyXestrous.di_met_time{nt}(na,:);
        xData(2,:) = b(2).XData+b(2).XOffset;
        xData(1,:) = b(1).XData+b(1).XOffset;


        
        plot(xData,aPlot,'o-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.25)
        clear aPlot xData
        
    end
end
legend(b,{'di/met';'est/pro'})
ylabel('Trial initiation latency (z-score)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])






%% Plot by time in session, z-score by animal, strict classification 
clear latencyXestrous 
Latency_Z = nan(size(Latency));
timeThresh = 3600;
for na = 1:numel(ids)
      idx = find(strcmp(Animal,ids{na})&TrialTime<=timeThresh&Dates>=sessionStart);
      Latency_Z(idx) = nanzscore(Latency(idx));
end

groups = {'D';'P';'E';'M'};
%groups = {'D';'D/P';'P';'P/E';'E';'E/M';'M';'M/D'};

[~,fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC);
fmap = fmap(round(linspace(20,size(fmap,1),numel(groups))),:);

timeBins = linspace(0,3600,2);
[~,~,TimeBins] = histcounts(TrialTime,timeBins);
for na = 1:numel(ids)
    for ntb = 1:numel(timeBins)-1
        
        for nb = 1:binNum
            for ng = 1:numel(groups)
                if numel(groups{ng})>1
                    gname = cat(2,groups{ng}(1), groups{ng}(end));
                else
                    gname = groups{ng};
                end
                eval(sprintf('latencyXestrous.%s{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&strcmp(Estrous_strict,''%s'')));',gname,groups{ng}));
            end
            latencyXestrous.est_pro_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous_strict,'P')|strcmp(Estrous_strict,'E'))));
            latencyXestrous.di_met_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&(strcmp(Estrous_strict,'M')|strcmp(Estrous_strict,'D'))));
            
            latencyXestrous.overall_time{ntb}(na,nb) = nanmean(Latency_Z(TimeBins==ntb&Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb));
            
        end
    end
end
clear thisPlot
figure();
clear bb
for nt = 1:numel(timeBins)-1
    for ng = 1:numel(groups)
        if numel(groups{ng})>1
            gname = cat(2,groups{ng}(1), groups{ng}(end));
        else
            gname = groups{ng};
        end
        thisPlot(ng,:) = nanmean(eval(sprintf('latencyXestrous.%s{nt}',gname)));
        thisError(ng,:) = nanstd(eval(sprintf('latencyXestrous.%s{nt}',gname)))./sqrt(sum(~isnan(eval(sprintf('latencyXestrous.%s{nt}(:,1)',gname)))));
    end
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    pause(0.1)
    for ng = 1:numel(groups)
        b(ng).EdgeColor = fmap(ng,:);
        b(ng).FaceColor = 'none';
        b(ng).LineWidth = 2;
      %  b(ng).FaceAlpha = .7;
    end
    for na = 1:numel(ids)
        clear aPlot xData
        for ng = 1:numel(groups)
            if numel(groups{ng})>1
                gname = cat(2,groups{ng}(1), groups{ng}(end));
            else
                gname = groups{ng};
            end
            aPlot(ng,:) = eval(sprintf('latencyXestrous.%s{nt}(na,:)',gname));
            xData(ng,:) = b(ng).XData+b(ng).XOffset;
            bb(ng) = plot(xData(ng,:),aPlot(ng,:),'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none');
        end
        
       p=plot(xData,aPlot,'-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
       uistack(bb,'top') 
        
    end
end
ylabel('Trial initiation latency (s)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])
legend(groups)

clear thisPlot
figure();
for nt = 1:numel(timeBins)-1
    
    thisPlot(2,:) = nanmean(latencyXestrous.est_pro_time{nt});
    thisPlot(1,:) = nanmean(latencyXestrous.di_met_time{nt});
   
    subplot(1,numel(timeBins)-1,nt); hold on
    b=bar(thisPlot');
    for ng = 1:2
        b(ng).EdgeColor = b(ng).FaceColor;
        b(ng).FaceColor = 'none';
    end
    for na = 1:numel(ids)
        aPlot(2,:) = latencyXestrous.est_pro_time{nt}(na,:);
        aPlot(1,:) = latencyXestrous.di_met_time{nt}(na,:);
        xData(2,:) = b(2).XData+b(2).XOffset;
        xData(1,:) = b(1).XData+b(1).XOffset;


        
        plot(xData,aPlot,'o-','Color',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'LineWidth',.25)
        clear aPlot xData
        
    end
end
legend(b,{'di/met';'est/pro'})
ylabel('Trial initiation latency (z-score)')
xlabel('Relative chosen value')
set(gca,'XTick',[1:binNum])





clear aPlot
f=figure(); hold on
f.Position = [306   463   167   252];
for ng = 1:numel(groups)
    thisdata = eval(sprintf('numTrials.%s',groups{ng}));
   thismean = nanmean(thisdata);
   thissem = nanstd(eval(sprintf('numTrials.%s',groups{ng})))./sqrt(sum(~isnan(eval(sprintf('numTrials.%s',groups{ng}))))); 
%    p = plot(ng,thismean,'o','Color',fmap(ng,:));
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
   plot([ng ng], [thismean-thissem, thismean+thissem], 'Color',fmap(ng,:),'LineWidth',1.5)
   plot([ng-.1 ng+.1], [thismean, thismean], 'Color',fmap(ng,:),'LineWidth',1.5)

   plot(repmat(ng+.25,1,numel(thisdata)),thisdata,'o','MarkerFaceColor',fmap(ng,:), 'MarkerEdgeColor','none','MarkerSize',5);
   aPlot(ng,:) = thisdata;
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Number of trials')
xlabel('Estrous stage')
set(gca,'FontSize',10)
clear p
p=plot(repmat([1:numel(groups)]'+.25, 1, numel(ids)),aPlot,'Color',[.7 .7 .7],'LineWidth',.5);



nback = 5;
for na = 1:numel(ids)
    thisSess = unique(Session(strcmp(Animal,ids{na})&Dates>=sessionStart));
    sessIdx = arrayfun(@(x) find(Session==x&strcmp(Animal,ids{na}),1,'first'),thisSess);
    thisSess = thisSess(strcmp(Estrous(sessIdx),'E')|strcmp(Estrous(sessIdx),'P'))';
    thisChoice = Choice(strcmp(Animal,ids{na})&~isnan(Latency_thresh)&boolean(sum(Session==thisSess,2)));
    thisChoice(thisChoice==-1) = NaN;
    thisStay = [NaN; thisChoice(1:end-1)==thisChoice(2:end)];
    thisPrev = [];
    thisRew = [];
    tempTrialCount = [];
    for ns = 1:numel(thisSess)
        temp   = Reward(strcmp(Animal,ids{na})&Session==thisSess(ns)&~isnan(Latency_thresh));
        thisPrev = cat(1,thisPrev,[NaN; temp(1:end-1)]);
        thisRew  = cat(1, thisRew,temp);
    end
   numTrials.EP(na) = mean(tempTrialCount);
   numTrials.EP_sem(na) = std(tempTrialCount)./sqrt(numel(tempTrialCount));
    
    choiceXoutcome.prevRew_EP(na) = nanmean(thisStay(thisPrev==1));
    choiceXoutcome.prevNRew_EP(na)= nanmean(thisStay(thisPrev==0));
    
    if isempty(thisSess)
        coeffs.EP(:,na) = nan(nback*2 + 1,1);
    else
    R = nan(size(thisChoice));
    UR = nan(size(thisChoice));
    
    %Reward
    
    R(thisChoice == 0 & thisRew == 1,1) = -1;
    R(thisChoice == 0 & thisRew == 0,1) = 0;
    UR(thisChoice == 0 & thisRew == 0,1) = -1;
    UR(thisChoice == 0 & thisRew == 1,1) = 0;
    R(thisChoice == 1 & thisRew == 1,1) = 1;
    R(thisChoice == 1 & thisRew == 0,1) = 0;
    UR(thisChoice == 1 & thisRew == 0,1) = 1;
    UR(thisChoice == 1 & thisRew == 1,1) = 0;
    
    resp = thisChoice;
    resp = resp(2:end);
    
    clear M
    for i = 1:nback
        M(:,i) = [zeros(i-1,1).*NaN; R(1:end-(i-1)-1)];
        M(:,nback+i) = [zeros(i-1,1).*NaN; UR(1:end-(i-1)-1)];
    end
    mdl  = glmfit(M,resp,'binomial');
    coeffs.EP(:,na) = mdl;
    end
    thisSess = min(Session(strcmp(Animal,ids{na}))):max(Session(strcmp(Animal,ids{na})));
    sessIdx = arrayfun(@(x) find(Session==x&strcmp(Animal,ids{na}),1,'first'),thisSess);
    thisSess = thisSess(strcmp(Estrous(sessIdx),'D')|strcmp(Estrous(sessIdx),'M'));
    thisChoice = Choice(strcmp(Animal,ids{na})&boolean(sum(Session==thisSess,2)));
    thisChoice(thisChoice==-1) = NaN;
    thisStay = [NaN; thisChoice(1:end-1)==thisChoice(2:end)];
    thisPrev = [];
    thisRew = [];
    tempTrialCount = [];
    for ns = 1:numel(thisSess)
        temp   = Reward(strcmp(Animal,ids{na})&Session==thisSess(ns));
        thisPrev = cat(1,thisPrev,[NaN; temp(1:end-1)]);
        thisRew  = cat(1, thisRew,temp);
        tempTrialCount(ns) = numel(temp); 
    end
   numTrials.DM(na) = mean(tempTrialCount);
   numTrials.DM_sem(na) = std(tempTrialCount)./sqrt(numel(tempTrialCount));
    
    
    choiceXoutcome.prevRew_DM(na) = nanmean(thisStay(thisPrev==1));
    choiceXoutcome.prevNRew_DM(na)= nanmean(thisStay(thisPrev==0));
    
    if isempty(thisSess)
        coeffs.DM(:,na) = nan(size(mdl));
    else
    R = nan(size(thisChoice));
    UR = nan(size(thisChoice));
    
    %Reward
    
    R(thisChoice == 0 & thisRew == 1,1) = -1;
    R(thisChoice == 0 & thisRew == 0,1) = 0;
    UR(thisChoice == 0 & thisRew == 0,1) = -1;
    UR(thisChoice == 0 & thisRew == 1,1) = 0;
    R(thisChoice == 1 & thisRew == 1,1) = 1;
    R(thisChoice == 1 & thisRew == 0,1) = 0;
    UR(thisChoice == 1 & thisRew == 0,1) = 1;
    UR(thisChoice == 1 & thisRew == 1,1) = 0;
    
    resp = thisChoice;
    resp = resp(2:end);
    
    clear M
    for i = 1:nback
        M(:,i) = [zeros(i-1,1).*NaN; R(1:end-(i-1)-1)];
        M(:,nback+i) = [zeros(i-1,1).*NaN; UR(1:end-(i-1)-1)];
    end
    mdl  = glmfit(M,resp,'binomial');
    coeffs.DM(:,na) = mdl;
    end
    
    
      
end



groups = {'EP';'DM'};
figure(); hold on
for ng = 1:numel(groups)
    mdl = eval(sprintf('coeffs.%s',groups{ng}));
    p(ng)=errorbar(1:nback, nanmean(mdl(2:nback+1,:),2),nanstd(mdl(2:nback+1,:),[],2)./sqrt(numel(ids)),'CapSize',0,'LineWidth',1.5);
    errorbar(1:nback, nanmean(mdl(nback+2:end,:),2),nanstd(mdl(nback+2:end,:),[],2)./sqrt(numel(ids)),'LineStyle','--','CapSize',0,'LineWidth',1.5,'Color',p(ng).Color)
    plot([0 nback+1],[0 0],'Color',[.7 .7 .7])
    % for na = 1:numel(ids)
    %     plot(1:nback, mdl(2:nback+1,na),'Color',[0 0 0 .5])
    %     plot(1:nback, mdl(nback+2:end,na),'--','Color',[0 0 0 .5])
    % end
    ylabel('Coefficient')
    xlabel('Trial back')
end
legend(p,groups);

f=figure('Units','inches','Position',[5,5,6 2.5]); hold on
for ng =1:numel(groups)
    subplot(1,numel(groups),ng), hold on
    thischoice.prevRew = eval(sprintf('choiceXoutcome.prevRew_%s',groups{ng}));
    thischoice.prevNRew = eval(sprintf('choiceXoutcome.prevNRew_%s',groups{ng}));
    b = bar([nanmean(thischoice.prevRew), nanmean(thischoice.prevNRew)]);
    b.EdgeColor = b.FaceColor;
    b.FaceColor = 'none';
    b.LineWidth = 2;
    errorbar(1:2,[nanmean(thischoice.prevRew), nanmean(thischoice.prevNRew)],[nanstd(thischoice.prevRew)./sqrt(sum(~isnan(thischoice.prevNRew))), nanstd(thischoice.prevNRew)./sqrt(sum(~isnan(thischoice.prevNRew)))],'Color',b.EdgeColor,'CapSize',0,'LineStyle','none')
    %  hold on, plot([choiceXoutcome.prevRew; choiceXoutcome.prevNRew],'Color',[0 0 0 .5])
    ylabel('Stay probability')
    if ng == 2
        xlabel('Previous outcome')
    end
    set(gca,'XTick',[1 2],'XTickLabel',{'Rew';'No rew'},'YLim',[0 1])
    title(groups{ng})
end
%%

groups = {'D';'P';'E';'M'};
% Extract random effects for each animal 
[randEff,randEffNames] = randomEffects(mdl);
for na = 1:numel(ids)
   thisIntercept = randEff(strcmp(randEffNames.Level,ids{na})&strcmp(randEffNames.Name,'(Intercept)'))+mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'(Intercept)')); 
   thisVal = randEff(strcmp(randEffNames.Level,ids{na})&strcmp(randEffNames.Name,'Value'))+mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'Value'));
   thisSess = randEffNames.Level(strcmp(randEffNames.Group,'Session:Animal') & contains(randEffNames.Level,ids{na}));
   thisSess = unique(cellfun(@(x) str2num(x(1:strfind(x,' '))), thisSess));
   thisIdx = find(strcmp(randEffNames.Group,'Session:Animal') & contains(randEffNames.Level,ids{na})&strcmp(randEffNames.Name,'Value'));
   thisIdx_inter = find(strcmp(randEffNames.Group,'Session:Animal') & contains(randEffNames.Level,ids{na})&strcmp(randEffNames.Name,'(Intercept)'));
   thisEstrous = arrayfun(@(x) Estrous(find(Session==x&strcmp(Animal,ids{na}),1,'first')), thisSess);
   for ng = 1:numel(groups)
       eval(sprintf('valueCoeff.%s(na) = thisVal+mean(randEff(thisIdx(strcmp(thisEstrous,groups{ng}))));',groups{ng}));
       eval(sprintf('valueCoeff.%s_all{na} = thisVal + randEff(thisIdx(strcmp(thisEstrous,groups{ng})));',groups{ng}));   
       eval(sprintf('valueCoeff.%s_intercept(na) = thisIntercept + mean(randEff(thisIdx_inter(strcmp(thisEstrous,groups{ng}))));',groups{ng}));
   end
end

clear thisVal
figure(); hold on
for ng = 1:numel(groups)
   thisVal(ng,:) = eval(sprintf('valueCoeff.%s',groups{ng}));
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Value coefficient')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5]) 
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end
clear thisVal
figure(); hold on
for ng = 1:numel(groups)
   thisVal(ng,:) = eval(sprintf('valueCoeff.%s_intercept',groups{ng}));
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Intercept')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5]) 
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end



%% Fit regression individually 
X = table(categorical(Session), Latency, (Estrous), Value, (Animal),TrialTime,Dates,'VariableNames',{'Session','Latency','Estrous','Value','Animal','TrialTime','Dates'});
X = X(strcmp(X.Estrous,'D')|strcmp(X.Estrous,'M')|strcmp(X.Estrous,'P')|strcmp(X.Estrous,'E'),:);
X.Estrous = categorical(X.Estrous);
X = X(X.Dates>=sessionStart,:);
X = X(X.TrialTime<=3600,:);

groups = {'D';'P';'E';'M'};
f = 'Latency ~ Value + (Value|Session)';

for na = 1:numel(ids)
    idx   = strcmp(X.Animal,(ids{na}));
    thisX = X(idx,:);
    mdl = fitglme(thisX,f,'DummyVarCoding','effect');

    
    % Extract random effects for each animal
    [randEff,randEffNames] = randomEffects(mdl);
   thisIntercept = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'(Intercept)')); 
   thisVal = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'Value')); 
  
   thisSess = unique(cellfun(@(x) str2num(x), randEffNames.Level(strcmp(randEffNames.Group,'Session'))));
   thisEstrous = arrayfun(@(x) Estrous(find(Session==x&strcmp(Animal,ids{na}),1,'first')), thisSess);
   thisIdx = find(strcmp(randEffNames.Name,'Value'));
   thisIdx_inter = find(strcmp(randEffNames.Name,'(Intercept)'));
   for ng = 1:numel(groups)
       eval(sprintf('valueCoeff.%s_indFit(na) = thisVal+mean(randEff(thisIdx(strcmp(thisEstrous,groups{ng}))));',groups{ng}));
       eval(sprintf('valueCoeff.%s_all_indFit{na} = thisVal + randEff(thisIdx(strcmp(thisEstrous,groups{ng})));',groups{ng}));
       eval(sprintf('valueCoeff.%s_intercept_indFit(na) = thisIntercept + mean(randEff(thisIdx_inter(strcmp(thisEstrous,groups{ng}))));',groups{ng}));
   end
end
clear thisVal
figure(); hold on
for ng = 1:numel(groups)
   thisVal(ng,:) = eval(sprintf('valueCoeff.%s_indFit',groups{ng}));
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
  p(ng)= errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Value coefficient')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5],'LineWidth',1) 
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end
for ng = 1:numel(groups)
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'Color',p(ng).Color,'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end


clear thisVal
figure(); hold on


for ng = 1:numel(groups)
   thisVal(ng,:) = eval(sprintf('valueCoeff.%s_intercept_indFit',groups{ng}));
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   p(ng)=errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Intercept')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5],'LineWidth',1) 
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end

for ng = 1:numel(groups)
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'Color',p(ng).Color,'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end


%% Fit regression individually 
X = table(categorical(Session), Latency, (Estrous), Value, (Animal),TrialTime,Dates,'VariableNames',{'Session','Latency','Estrous','Value','Animal','TrialTime','Dates'});
X = X(strcmp(X.Estrous,'D')|strcmp(X.Estrous,'M')|strcmp(X.Estrous,'P')|strcmp(X.Estrous,'E'),:);

X.Estrous = categorical(X.Estrous);
X = X(X.Dates>=sessionStart,:);
%X = X(X.TrialTime<=3600,:);

groups = {'D';'D/P';'P';'P/E';'E';'E/M';'M';'M/D'};
groups = {'D';'P';'E';'M'};
f = 'Latency ~ Value + (Value|Session)';

for na = 1:numel(ids)
    idx   = strcmp(X.Animal,(ids{na}));
    thisX = X(idx,:);
    mdl = fitlme(thisX,f,'DummyVarCoding','effect');
    
    % Extract random effects for each animal
    [randEff,randEffNames] = randomEffects(mdl);
   thisIntercept = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'(Intercept)')); 
   thisVal = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'Value')); 
  
   thisSess = unique(cellfun(@(x) str2num(x), randEffNames.Level(strcmp(randEffNames.Group,'Session'))));
   thisEstrous = arrayfun(@(x) Estrous(find(Session==x&strcmp(Animal,ids{na}),1,'first')), thisSess);
   thisIdx = find(strcmp(randEffNames.Name,'Value'));
   thisIdx_inter = find(strcmp(randEffNames.Name,'(Intercept)'));
   for ng = 1:numel(groups)
       if numel(groups{ng})>1
           gname = cat(2,groups{ng}(1),groups{ng}(end));
       else
           gname = groups{ng};
       end
       eval(sprintf('valueCoeff.%s_indFit(na) = thisVal+mean(randEff(thisIdx(strcmp(thisEstrous,groups{ng}))));',gname));
       eval(sprintf('valueCoeff.%s_all_indFit{na} = thisVal + randEff(thisIdx(strcmp(thisEstrous,groups{ng})));',gname));
       eval(sprintf('valueCoeff.%s_intercept_indFit(na) = thisIntercept + mean(randEff(thisIdx_inter(strcmp(thisEstrous,groups{ng}))));',gname));
   end
end
clear thisVal
figure(); hold on
for ng = 1:numel(groups)
    if numel(groups{ng})>1
        gname = cat(2,groups{ng}(1),groups{ng}(end));
    else
        gname = groups{ng};
    end
    thisVal(ng,:) = eval(sprintf('valueCoeff.%s_indFit',gname));
    %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
    p(ng)= errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1.5);
    %    p.MarkerFaceColor = p.Color;
    %    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Value coefficient')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5])
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end
for ng = 1:numel(groups)
    %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
    errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'Color',p(ng).Color,'LineStyle','none','LineWidth',1.5);
    %    p.MarkerFaceColor = p.Color;
    %    p.MarkerEdgeColor = 'none';
end


clear thisVal
figure(); hold on


for ng = 1:numel(groups)
      if numel(groups{ng})>1
        gname = cat(2,groups{ng}(1),groups{ng}(end));
    else
        gname = groups{ng};
    end
 
   thisVal(ng,:) = eval(sprintf('valueCoeff.%s_intercept_indFit',gname));
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   p(ng)=errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
set(gca,'XLim', [.5 ng+.5],'XTick',1:numel(groups),'XTickLabel',groups)
ylabel('Intercept')
xlabel('Estrous stage')

for na = 1:numel(ids)
    plot(1:ng,thisVal(:,na),'Color',[.7 .7 .7 .5],'LineWidth',1) 
    scatter(1:ng,thisVal(:,na), 40,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerFaceAlpha',.8)
end

for ng = 1:numel(groups)
   %p=plot(ng, nanmean(thisVal(ng,:)), 'o');
   errorbar([ng],nanmean(thisVal(ng,:)), nanstd(thisVal(ng,:))./sqrt(sum(~isnan(thisVal(ng,:)))),'Color',p(ng).Color,'LineStyle','none','LineWidth',1.5);
%    p.MarkerFaceColor = p.Color;
%    p.MarkerEdgeColor = 'none';
end
%% trashcan
%% Plot individuals
% nrows = ceil(numel(ids)/4);
% groups= {'D';'P';'E';'M'};
% 
% figure()
% for na = 1:numel(ids)
%     subplot(nrows,4,na); hold on
%     clear thisMean thisStd
%     for ng = 1:numel(groups)
%         for nb = 1:binNum
%             idx = find(Dates>=sessionStart&strcmp(Animal,ids{na})&Value_ptile==nb&TrialTime<=3600&strcmp(Estrous,groups{ng}));
%             thisMean(ng,nb) = nanmean(Latency(idx));
%             thisStd(ng,nb) = nanstd(Latency(idx))./sqrt(numel(idx)); 
%         end
%     end
%     errorbar(repmat(1:binNum,4,1)',thisMean',thisStd','CapSize',0)
%     set(gca,'XLim',[0.5 binNum+.5])
%     if na == 1
%         legend(groups)
%     end
% end
% 
% stats.individual = table(nan(numel(ids),1),nan(numel(ids),1),nan(numel(ids),1),nan(numel(ids),1),nan(numel(ids),1),nan(numel(ids),1),nan(numel(ids),1),'VariableNames',{'Value','Estrous','ValueXEstrous','Weight','ValueXWeight','EstrousXWeight','ValueXEstrousXWeight'});
% for na = 1:numel(ids)
%     idx = find(strcmp(Animal,ids{na}) & TrialTime<3600&Dates>=sessionStart);
%     X = table(Latency(idx),Session(idx),categorical(Value_ptile(idx)),categorical(Estrous(idx)),zscore(Weight(idx)),'VariableNames',{'Latency';'Session';'Value';'Estrous';'Weight'});
%     f = 'Latency ~ Value*Estrous*Weight + (1|Session)';
%     try
%         mdl = fitlme(X,f,'DummyVarCoding','effect');
%         a = anova(mdl,'DFMethod','Satterthwaite');
%         stats.individual.Value(na) = a.pValue(strcmp(a.Term,'Value'));
%         stats.individual.Estrous(na) = a.pValue(strcmp(a.Term,'Estrous'));
%         stats.individual.ValueXEstrous(na) = a.pValue(strcmp(a.Term,'Value:Estrous'));
%         stats.individual.Weight(na) = a.pValue(strcmp(a.Term,'Weight'));
%         stats.individual.ValueXWeight(na) = a.pValue(strcmp(a.Term,'Value:Weight'));
%         stats.individual.EstrousXWeight(na) = a.pValue(strcmp(a.Term,'Estrous:Weight'));
%         stats.individual.ValueXEstrousXWeight(na) = a.pValue(strcmp(a.Term,'Value:Estrous:Weight'));
%     catch
%         try
%             f = 'Latency ~ Value*Estrous + (1|Session)';
%             
%             mdl = fitlme(X,f,'DummyVarCoding','effect');
%             a = anova(mdl,'DFMethod','Satterthwaite');
%             stats.individual.Value(na) = a.pValue(strcmp(a.Term,'Value'));
%             stats.individual.Estrous(na) = a.pValue(strcmp(a.Term,'Estrous'));
%             stats.individual.ValueXEstrous(na) = a.pValue(strcmp(a.Term,'Value:Estrous'));
%         catch
%         end
%     end
% end
