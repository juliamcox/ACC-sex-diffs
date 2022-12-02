function humanBandit_valueXage_choice(choice_f,choice_m,valType,binNum)
age = choice_f.age;
age = cat(2,age,choice_m.age);

ageBins = prctile(age,linspace(0,100,3));
[~,~,age_f] = histcounts(choice_f.age,ageBins);
[~,~,age_m] = histcounts(choice_m.age,ageBins);

for nv = 1:numel(valType)
        thisChoice = eval(sprintf('choice_f.probRight_%s;',valType{nv}));
        for nb = 1:numel(ageBins)
            thisIDs = find(age_f==nb);
            for na = 1:numel(thisIDs)
                choicePlot{nb}(na,:) = thisChoice(thisIDs(na),:);
            end
        end
end

mu_f = cellfun(@(x) nanmean(x), choicePlot,'UniformOutput',false);
sem_f = cellfun(@(x) nanstd(x)./sqrt(size(x,1)),choicePlot,'UniformOutput',false);
for nv = 1:numel(valType)
        thisChoice = eval(sprintf('choice_m.probRight_%s;',valType{nv}));
        for nb = 1:numel(ageBins)
            thisIDs = find(age_m==nb);
            for na = 1:numel(thisIDs)
                choicePlot_m{nb}(na,:) = thisChoice(thisIDs(na),:);
            end
        end
end

mu_m = cellfun(@(x) nanmean(x), choicePlot_m,'UniformOutput',false);
sem_m = cellfun(@(x) nanstd(x)./sqrt(size(x,1)),choicePlot_m,'UniformOutput',false);
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));


figure(), hold on
for nb = 1:numel(ageBins)-1
    subplot(1,2,nb); hold on
    errorbar([1:binNum],mu_f{nb},sem_f{nb},'Color',plotParams.femaleC,'LineWidth',1.5,'CapSize',0);
    
    axis square
    set(gca,'XLim',[.5 binNum+.5])
    
    errorbar([1:binNum],mu_m{nb},sem_m{nb},'Color',plotParams.maleC,'LineWidth',1.5,'CapSize',0);
    
    axis square
    set(gca,'YLim', [0 1])
    ylabel('P(choice=r)')
    xlabel(valType{1})
end
set(gca,'XLim',[.5 binNum+.5])