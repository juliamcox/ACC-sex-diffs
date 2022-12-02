function postHoc = humanBandit_valueXage(latency_f,latency_m,valType,latencyType,binNum,zscoreFlag)

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

age = latency_f.age;
age = cat(2,age,latency_m.age);

ageBins = prctile(age,linspace(0,100,3));
[~,~,age_f] = histcounts(latency_f.age,ageBins);
[~,~,age_m] = histcounts(latency_m.age,ageBins);

for nv = 1:numel(valType)
    for nl = 1:numel(latencyType)
        thisLatency = eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv}));
        for nb = 1:numel(ageBins)
            thisIDs = find(age_f==nb);
            for na = 1:numel(thisIDs)
                latencyPlot{nb}(na,:) = thisLatency(thisIDs(na),:);
            end
        end
    end
end

mu_f = cellfun(@(x) nanmean(x), latencyPlot,'UniformOutput',false);
sem_f = cellfun(@(x) nanstd(x)./sqrt(size(x,1)),latencyPlot,'UniformOutput',false);
for nv = 1:numel(valType)
    for nl = 1:numel(latencyType)
        thisLatency = eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv}));
        for nb = 1:numel(ageBins)
            thisIDs = find(age_m==nb);
            for na = 1:numel(thisIDs)
                latencyPlot_m{nb}(na,:) = thisLatency(thisIDs(na),:);
            end
        end
    end
end

mu_m = cellfun(@(x) nanmean(x), latencyPlot_m,'UniformOutput',false);
sem_m = cellfun(@(x) nanstd(x)./sqrt(size(x,1)),latencyPlot_m,'UniformOutput',false);

for nb = 1:numel(ageBins)-1
    for nbb = 1:size(latencyPlot_m{nb},2)
        [postHoc.p{nb}(nbb), ~,  postHoc.stats{nb}(nbb)] = ranksum(latencyPlot_m{nb}(:,nbb),latencyPlot{nb}(:,nbb));
    end
end


figure(), hold on
for nb = 1:numel(ageBins)-1
    subplot(1,numel(ageBins)-1,nb); hold on
    errorbar([1:binNum],mu_f{nb},sem_f{nb},'Color',plotParams.femaleC,'LineWidth',1.5,'CapSize',0);
    if zscoreFlag
        set(gca,'YLim', [-0.05 0.09])
        ylabel('Trial initiation latency (Z-score)')
    else
        set(gca,'YLim', [200 650])
        ylabel('Trial initiation latency (s)')
    end
    axis square
    set(gca,'XLim',[.5 binNum+.5])
    
    errorbar([1:binNum],mu_m{nb},sem_m{nb},'Color',plotParams.maleC,'LineWidth',1.5,'CapSize',0);
    xlabel(valType{1})
end

