function latencyXweight(behaviorTable,binNum,ids_f_all,ids_m_all)


plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

aids = unique(behaviorTable.aID);

behaviorTable.trialInit(behaviorTable.trialInit==0,:) = 0.01; 
behaviorTable.weight_zscore = nanzscore(behaviorTable.weight); 

% find average weight per subject
behaviorTable.meanWeight = zeros(size(behaviorTable,1),1);
for na = 1:numel(aids)
   idx = find(strcmp(behaviorTable.aID, aids{na}));
   sessionIdx = unique(behaviorTable.session(idx));
   sessionIdx(isnan(sessionIdx)) = [];
   meanWeight = arrayfun(@(x) behaviorTable.weight(find(behaviorTable.session==x&strcmp(behaviorTable.aID, aids{na}),1,'first')),sessionIdx);
   behaviorTable.meanWeight(idx) = nanmean(meanWeight); 
end
behaviorTable.meanWeight_zscore = nanzscore(behaviorTable.meanWeight); 

%% Plot by weight terciles 
weightBinNum = 3; 

for na = 1:numel(ids_f_all)
   temp = behaviorTable(strcmp(behaviorTable.aID,ids_f_all{na})&behaviorTable.laserSession==0,:);
   thisSess = unique(temp.session);
   thisSess(isnan(thisSess)) = [];
   thisWeights = arrayfun(@(x) temp.weight(find(temp.session==x,1,'first')), thisSess); 
   thisWeightBins = prctile(thisWeights,linspace(0,100,weightBinNum+1));
  
   if sum(diff(thisWeightBins)==0)==0
       [~,~,bins] = histcounts(thisWeights,thisWeightBins);
       for nb = 1:weightBinNum
           sess = thisSess(bins==nb);
           for nbb = 1:binNum
               LatencyXvalueXweight{nb}(na,nbb) = nanmean(temp.trialInit_thresh(boolean(sum(temp.session==sess',2))&temp.qChosenDiff_quant==nbb),1);
           end
       end
   else
       for nb = 1:weightBinNum
           LatencyXvalueXweight{nb}(na,:) = nan(1,binNum);
       end
   end
end
for na = 1:numel(ids_m_all)
   temp = behaviorTable(strcmp(behaviorTable.aID,ids_m_all{na})&behaviorTable.laserSession==0,:);
   thisSess = unique(temp.session); 
   thisSess(isnan(thisSess)) = [];
   thisWeights = arrayfun(@(x) temp.weight(find(temp.session==x,1,'first')), thisSess); 
   thisWeightBins = prctile(thisWeights,linspace(0,100,weightBinNum+1));
   if sum(diff(thisWeightBins)==0)==0
       [~,~,bins] = histcounts(thisWeights,thisWeightBins);
       for nb = 1:weightBinNum
           sess = thisSess(bins==nb);
           for nbb = 1:binNum
               LatencyXvalueXweight_m{nb}(na,nbb) = nanmean(temp.trialInit_thresh(boolean(sum(temp.session==sess',2))&temp.qChosenDiff_quant==nbb),1);
           end
       end
   else
       for nb = 1:weightBinNum
           LatencyXvalueXweight_m{nb}(na,:) = nan(1,binNum);
       end
   end
end

[mmap, fmap] = maleFemaleColormap(plotParams.maleC,plotParams.femaleC,50);
figure('Position',[681 832 560 145]);
subplot(1,2,1); hold on
for nb = 1:weightBinNum
    
    p=shadedErrorBar(1:binNum,nanmean(LatencyXvalueXweight{nb}),nanstd(LatencyXvalueXweight{nb})./sqrt(sum(~isnan(LatencyXvalueXweight{nb}(:,1)),1)),'lineprops',{ 'Color',fmap((nb*10)+10,:)});
    %errorbar(1:binNum,nanmean(LatencyXvalueXweight_m{nb}),nanstd(LatencyXvalueXweight_m{nb})./sqrt(sum(~isnan(LatencyXvalueXweight_m{nb}(:,1)),1)), 'Color',plotParams.maleC,'CapSize',0)
    set(gca,'XLim',[0.5 binNum+.5],'YLim',[0 10])
    
    pp(nb) = p.mainLine;
end
legend(pp,{num2str([1:weightBinNum]')});

subplot(1,2,2); hold on
for nb = 1:weightBinNum
    
    p=shadedErrorBar(1:binNum,nanmean(LatencyXvalueXweight_m{nb}),nanstd(LatencyXvalueXweight_m{nb})./sqrt(sum(~isnan(LatencyXvalueXweight_m{nb}(:,1)),1)),'lineprops',{ 'Color',mmap((nb*10)+10,:)});
    %errorbar(1:binNum,nanmean(LatencyXvalueXweight_m{nb}),nanstd(LatencyXvalueXweight_m{nb})./sqrt(sum(~isnan(LatencyXvalueXweight_m{nb}(:,1)),1)), 'Color',plotParams.maleC,'CapSize',0)
    set(gca,'XLim',[0.5 binNum+.5],'YLim',[0 10])
    pp(nb) = p.mainLine;
end
legend(pp,{num2str([1:weightBinNum]')});





%% Plot by weight bin 
Weights = [];
Animal = [];
LatencyXvalue=[];
for na = 1:numel(ids_m_all)
    clear thisSess
    thisSess = unique(behaviorTable.session(behaviorTable.laserSession==0&strcmp(behaviorTable.aID,ids_m_all{na})));
    %thisSess(isnan(thisSess)) = [];
    temp = behaviorTable(strcmp(behaviorTable.aID,ids_m_all{na}),:);
    Weights = cat(1,Weights,arrayfun(@(x) temp.weight(find(temp.session==x,1,'first')),thisSess));
    Animal  = cat(1,Animal,repmat(ids_m_all(na),numel(thisSess),1));
    for ns = thisSess'
        for nb = 1:binNum
            thisVal(nb) = nanmean(temp.trialInit_thresh(temp.session==ns&temp.qChosenDiff_quant==nb));
        end
        LatencyXvalue = cat(1,LatencyXvalue,thisVal);
    end
end
for na = 1:numel(ids_f_all)
    clear thisSess
    thisSess = unique(behaviorTable.session(behaviorTable.laserSession==0&strcmp(behaviorTable.aID,ids_f_all{na})));
    temp = behaviorTable(strcmp(behaviorTable.aID,ids_f_all{na}),:);
    Weights = cat(1,Weights,arrayfun(@(x) temp.weight(find(temp.session==x,1,'first')),thisSess));
    Animal  = cat(1,Animal,repmat(ids_f_all(na),numel(thisSess),1));
    for ns = thisSess'
        for nb = 1:binNum
            thisVal(nb) = nanmean(temp.trialInit_thresh(temp.session==ns&temp.qChosenDiff_quant==nb));
        end
        LatencyXvalue = cat(1,LatencyXvalue,thisVal);
    end
end


weightBins = [14, 20,25,max(Weights)];
[~,~,binWeight] = histcounts(Weights,weightBins);

latencyXvalueXweight_m = cell(1,max(binWeight));
latencyXvalueXweight_f = cell(1,max(binWeight));

countXvalueXweight_m = cell(1,max(binWeight));
countXvalueXweight_f = cell(1,max(binWeight));

for nb = 1:numel(latencyXvalueXweight_m)
    latencyXvalueXweight_m{nb} = nan(numel(ids_m_all),binNum);
    latencyXvalueXweight_f{nb} = nan(numel(ids_f_all),binNum);
    countXvalueXweight_f{nb} = nan(numel(ids_f_all),1);
    countXvalueXweight_m{nb} = nan(numel(ids_m_all),1);
end

for na = 1:numel(ids_m_all)
   thisLatencyXval = LatencyXvalue(strcmp(Animal,ids_m_all{na}),:);
   thisBins = binWeight(strcmp(Animal,ids_m_all{na}));
   for nb = 1:numel(latencyXvalueXweight_m)
       latencyXvalueXweight_m{nb}(na,:) = nanmean(thisLatencyXval(thisBins==nb,:),1);
       countXvalueXweight_m{nb}(na,:) = nansum(thisBins==nb);
   end
end

for na = 1:numel(ids_f_all)
    thisLatencyXval = LatencyXvalue(strcmp(Animal,ids_f_all{na}),:);
    thisBins = binWeight(strcmp(Animal,ids_f_all{na}));
    for nb = 1:numel(latencyXvalueXweight_f)
        latencyXvalueXweight_f{nb}(na,:) = nanmean(thisLatencyXval(thisBins==nb,:),1);
        countXvalueXweight_f{nb}(na,:) = nansum(thisBins==nb);
    end
end

figure('Position',[681 832 560 145]);
for nb = 1:numel(latencyXvalueXweight_f)
    subplot(1,numel(latencyXvalueXweight_f),nb); hold on
%     errorbar(1:binNum,nanmean(latencyXvalueXweight_f{nb}), nanstd(latencyXvalueXweight_f{nb})./sqrt(sum(~isnan(latencyXvalueXweight_f{nb}(:,1)))),'Color',plotParams.femaleC,'CapSize',0,'LineWidth',1);
%     errorbar(1:binNum,nanmean(latencyXvalueXweight_m{nb}), nanstd(latencyXvalueXweight_m{nb})./sqrt(sum(~isnan(latencyXvalueXweight_m{nb}(:,1)))),'Color',plotParams.maleC,'CapSize',0,'LineWidth',1);
%     plot(repmat(1:binNum,size(latencyXvalueXweight_f{nb},1),1)',latencyXvalueXweight_f{nb}','Color',[plotParams.femaleC .5], 'LineWidth',.25) 
%     plot(repmat(1:binNum,size(latencyXvalueXweight_m{nb},1),1)',latencyXvalueXweight_m{nb}','Color',[plotParams.maleC .5], 'LineWidth',.25) 

    shadedErrorBar(1:binNum,nanmean(latencyXvalueXweight_f{nb}), nanstd(latencyXvalueXweight_f{nb})./sqrt(sum(~isnan(latencyXvalueXweight_f{nb}(:,1)))),'lineprops',{ 'Color',plotParams.femaleC});
    shadedErrorBar(1:binNum,nanmean(latencyXvalueXweight_m{nb}), nanstd(latencyXvalueXweight_m{nb})./sqrt(sum(~isnan(latencyXvalueXweight_m{nb}(:,1)))),'lineprops',{ 'Color',plotParams.maleC});

    set(gca,'YLim',[0 10],'XLim',[.5 binNum+.5])
    xlabel('Relative chosen value')
    ylabel('Trial initiation latency (s)')
   
    title(sprintf('%s - %s g', num2str(weightBins(nb)),num2str(weightBins(nb+1))))
end
