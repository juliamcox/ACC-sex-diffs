function stats=latencyQuant_optoXtime(cohort,behaviorTable,valType_plot,binNum,laserType,zscoreFlag,latencyType)


% Find opto mice 
ids_m = generateAnimalList(sprintf('%s_male',cohort));
ids_f = generateAnimalList(sprintf('%s_female',cohort));
ids_m_yfp = generateAnimalList(sprintf('%s_yfp_male',cohort));
ids_f_yfp = generateAnimalList(sprintf('%s_yfp_female',cohort));

behaviorTable = behaviorTable(behaviorTable.laserSession==1,:); 


thisValue_ptile = eval(sprintf('behaviorTable.%s_quant',valType_plot));
if contains(valType_plot,'thresh')
    thisValue = eval(sprintf('behaviorTable.%s',valType_plot(1:strfind(valType_plot,'_')-1)));
else
    thisValue = eval(sprintf('behaviorTable.%s',valType_plot));
end
behaviorTable.thisValue = thisValue;
behaviorTable.thisValue_ptile = thisValue_ptile;

if zscoreFlag
    Latency = eval(sprintf('behaviorTable.%s_zscore',latencyType));
else
    Latency = eval(sprintf('behaviorTable.%s',latencyType));
end
    
Latency(Latency<.01) = 0.01; % correct for different temporal resolutions

behaviorTable.Latency = Latency;
behaviorTable(isnan(behaviorTable.Latency),:) = [];

%% 
ntimeBins = 2;
rangeStart = 0:floor(7200/ntimeBins):7200;
if rangeStart(end) ~= 7200
    rangeStart(end+1) = 7201;
end

for nt = 1:ntimeBins
   X = behaviorTable(behaviorTable.opsin==1&behaviorTable.female==1&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)&(behaviorTable.laser==0|behaviorTable.laser==1),:);
   X.laser = categorical(X.laser);
   X.value = categorical(X.qChosenDiff_quant); 
   f = 'trialInit~laser*value + (1|aID)';
   mdl = fitlme(X,f,'DummyVarCoding','effect');
   stats_female{nt} = anova(mdl,'DFMethod','Satterthwaite');
end

for nt = 1:ntimeBins
   X = behaviorTable(behaviorTable.opsin==1&behaviorTable.female==0&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)&(behaviorTable.laser==0|behaviorTable.laser==1),:);
   X.laser = categorical(X.laser);
   X.value = categorical(X.qChosenDiff_quant); 
   f = 'trialInit~laser*value + (1|aID)';
   mdl = fitlme(X,f,'DummyVarCoding','effects');
   stats_male{nt} = anova(mdl,'DFMethod','Satterthwaite');
end

X = behaviorTable(behaviorTable.opsin==1&(behaviorTable.laser==0|behaviorTable.laser==1),:);
X.female = categorical(X.female);
X.value = categorical(X.qChosenDiff_quant);
X.trial = nanzscore(X.trial); 
X.laser = categorical(X.laser); 
X.timeBin = zeros(size(X,1),1);
for nt = 1:ntimeBins
    X.timeBin(X.trialTime>=rangeStart(nt)&X.trialTime<rangeStart(nt+1))=nt;
end
X(X.timeBin==0,:) = [];

X.timeBin = categorical(X.timeBin); 
f = 'trialInit~laser*value*trial*female + (1+value*trial|aID)';
mdl = fitlme(X,f,'DummyVarCoding','effects');

%% Extract animal means

groups = {'f';'m';'f_yfp';'m_yfp'};

for nt = 1:(ntimeBins)
for ng = 1:numel(groups)
    thisIDs = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIDs)
        for nb = 1:binNum
            eval(sprintf('X_%s{nt}(na,nb) = nanmean(behaviorTable.Latency(behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==0&behaviorTable.trial~=1&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)));',groups{ng}))
            eval(sprintf('X_laser_%s{nt}(na,nb) = nanmean(behaviorTable.Latency(behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==laserType&behaviorTable.trial~=1&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)));',groups{ng}))
            
            eval(sprintf('X_%s_count{nt}(na,nb) = nansum((behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==0&behaviorTable.trial~=1&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)));',groups{ng}))
            eval(sprintf('X_laser_%s_count{nt}(na,nb) = nansum((behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==laserType&behaviorTable.trial~=1&behaviorTable.trialTime>=rangeStart(nt)&behaviorTable.trialTime<rangeStart(nt+1)));',groups{ng}))
        end
    end
end


%% Female laser v no laser by value
for nb = 1:binNum
        whichTest{1}(nb) = 1;
        [stats{nt}.pvals_posthoc(nb,1),~,stats{nt}.posthoc{nb,1}] = signrank(X_f{nt}(:,nb), X_laser_f{nt}(:,nb));

end

% Male laser v no laser by value
for nb = 1:binNum
        whichTest{2}(nb) = 1;
        [stats{nt}.pvals_posthoc(nb,2),~,stats{nt}.posthoc{nb,2}] = signrank(X_m{nt}(:,nb), X_laser_m{nt}(:,nb));

end



end
%% Mixed-effects regession of  difference between laser and no laser
for nt = 1:ntimeBins
value = [];
LatencyDiff = [];
sex = [];
opsin = [];
animal = [];
for nb = 1:binNum
    LatencyDiff = cat(1,LatencyDiff,X_laser_f{nt}(:,nb)-X_f{nt}(:,nb),...
    X_laser_m{nt}(:,nb)-X_m{nt}(:,nb),...
    X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb),...
    X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb));

    value = cat(1,value,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))).*nb,...
    ones(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))).*nb,...
    ones(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))).*nb,...
    ones(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))).*nb);
    sex = cat(1,sex,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))),...
    zeros(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))),...
    ones(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))),...
    zeros(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))));
    opsin = cat(1,opsin,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))),...
    ones(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))),...
    zeros(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))),...
    zeros(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))));
    
    animal = cat(1,animal,[1:size(X_f{nt},1)]', [size(X_f{nt},1)+1:size(X_f{nt},1)+size(X_m{nt},1)]',[size(X_f{nt},1)+size(X_m{nt})+1:size(X_f{nt},1)+size(X_m{nt},1)+size(X_f_yfp{nt},1)]',...
        [size(X_f{nt},1)+size(X_m{nt})+size(X_f_yfp{nt})+1:size(X_f{nt},1)+size(X_m{nt},1)+size(X_f_yfp{nt},1)+size(X_m_yfp{nt},1)]');
    

end

T = table(LatencyDiff,'VariableNames',{'Latency'});
T.Sex = categorical(sex);
T.Opsin = categorical(opsin);
T.Value = categorical(value); 
T.Animal = categorical(animal);

f = 'Latency ~ Opsin*Value*Sex + (1|Animal)';
glme_diff = fitlme(T,f,'DummyVarCoding','effects');
stats{nt}.glme_diff_coeff = dataset2cell((glme_diff.Coefficients));
stats{nt}.glme_diff=glme_diff;
stats{nt}.glme_diff_anova = dataset2cell(anova(glme_diff,'DFMethod','Satterthwaite'));

T = behaviorTable;
T = T(T.laser==0|T.laser==1,:);
T = T(~isnan(T.trialInit_thresh),:);
T.trialInit_thresh(T.trialInit_thresh==0) = 0.01;

T.latency = log(T.trialInit_thresh);
T.female = categorical(T.female);
T.qChosenDiff_quant = categorical(T.qChosenDiff_quant); 
T.opsin = categorical(T.opsin);
T.laser = categorical(T.laser);



end

% stats
    
value = [];
LatencyDiff = [];
sex = [];
opsin = [];
animal = [];
sessiontime = [];
for nt = 1:ntimeBins

for nb = 1:binNum
    LatencyDiff = cat(1,LatencyDiff,X_laser_f{nt}(:,nb)-X_f{nt}(:,nb),...
    X_laser_m{nt}(:,nb)-X_m{nt}(:,nb),...
    X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb),...
    X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb));

    value = cat(1,value,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))).*nb,...
    ones(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))).*nb,...
    ones(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))).*nb,...
    ones(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))).*nb);
    sex = cat(1,sex,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))),...
    zeros(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))),...
    ones(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))),...
    zeros(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))));
    opsin = cat(1,opsin,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))),...
    ones(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))),...
    zeros(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))),...
    zeros(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))));

    sessiontime = cat(1,sessiontime,ones(size(X_laser_f{nt}(:,nb)-X_f{nt}(:,nb))).*nt,...
    ones(size(X_laser_m{nt}(:,nb)-X_m{nt}(:,nb))).*nt,...
    ones(size(X_laser_f_yfp{nt}(:,nb)-X_f_yfp{nt}(:,nb))).*nt,...
    ones(size(X_laser_m_yfp{nt}(:,nb)-X_m_yfp{nt}(:,nb))).*nt);

    animal = cat(1,animal,[1:size(X_f{nt},1)]', [size(X_f{nt},1)+1:size(X_f{nt},1)+size(X_m{nt},1)]',[size(X_f{nt},1)+size(X_m{nt})+1:size(X_f{nt},1)+size(X_m{nt},1)+size(X_f_yfp{nt},1)]',...
        [size(X_f{nt},1)+size(X_m{nt})+size(X_f_yfp{nt})+1:size(X_f{nt},1)+size(X_m{nt},1)+size(X_f_yfp{nt},1)+size(X_m_yfp{nt},1)]');
    

end
end
T = table(LatencyDiff,'VariableNames',{'Latency'});
T.Sex = categorical(sex);
T.Opsin = categorical(opsin);
T.Value = categorical(value); 
T.Animal = categorical(animal);
T.time = categorical(sessiontime); 

f = 'Latency ~ Opsin*Value*Sex*time + (1+Value|Animal)';
glme_diff = fitlme(T,f,'DummyVarCoding','effects');
stats{nt+1}.glme_diff_coeff = dataset2cell((glme_diff.Coefficients));
stats{nt+1}.glme_diff=glme_diff;
stats{nt+1}.glme_diff_anova = dataset2cell(anova(glme_diff,'DFMethod','Satterthwaite'));




%% Plot 
    
    
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;

   
  


figure();
% Plot female nphr outcome laser v no laser
for nt = 1:(ntimeBins)
    subplot(2,ntimeBins,1+(nt-1)); hold on
    for nb = 1:binNum
        xaxis = cat(2,repmat(nb,1,size(X_f{nt},1))'-.2,repmat(nb,1,size(X_f{nt},1))'+.2);
        plot(xaxis',cat(2,X_f{nt}(:,nb),X_laser_f{nt}(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
    end
    for nb = 1:binNum
        xaxis = repmat(nb,1,size(X_f{nt},1))-.2;%+datasample(0:.01:.25,size(X_f,1));
        scatter(xaxis, X_f{nt}(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
    end
    for nb = 1:binNum
        xaxis = repmat(nb,1,size(X_f{nt},1))+.2;%-datasample(0:.01:.25,size(X_f,1));
        scatter(xaxis, X_laser_f{nt}(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
    end
    
    p=errorbar(nanmean(X_laser_f{nt}),nanstd(X_laser_f{nt})./sqrt(size(X_laser_f{nt},1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0);
    p(2)=errorbar(1:binNum,nanmean(X_f{nt}),nanstd(X_f{nt})./sqrt(size(X_f{nt},1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
    legend(p,{'laser';'no laser'})
    
    subplot(2,ntimeBins,ntimeBins+(nt)); hold on
    for nb = 1:binNum
        xaxis = cat(2,repmat(nb,1,size(X_m{nt},1))'-.2,repmat(nb,1,size(X_m{nt},1))'+.2);
        plot(xaxis',cat(2,X_m{nt}(:,nb),X_laser_m{nt}(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
    end
    for nb = 1:binNum
        xaxis = repmat(nb,1,size(X_m{nt},1))-.2;%+datasample(0:.01:.25,size(X_m,1));
        scatter(xaxis, X_m{nt}(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
    end
    for nb = 1:binNum
        xaxis = repmat(nb,1,size(X_m{nt},1))+.2;%-datasample(0:.01:.25,size(X_m,1));
        scatter(xaxis, X_laser_m{nt}(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
    end
    p = errorbar(nanmean(X_m{nt}),nanstd(X_m{nt})./sqrt(size(X_m{nt},1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
    p(2)=errorbar(nanmean(X_laser_m{nt}),nanstd(X_laser_m{nt})./sqrt(size(X_laser_m{nt},1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0);
    legend(p,{'laser';'no laser'})
        
    for np = 1:ntimeBins*2
        subplot(2,ntimeBins,np)
        axis square
        ylabel('Trial initiation latency (sec)')
        set(gca,'YLim',[0,20],'XLim', [.5 max(unique(behaviorTable.thisValue_ptile))+.5])
        if laserType==3
            set(gca,'YLim',[0 20],'XLim', [.5 max(unique(behaviorTable.thisValue_ptile))+.5])
        end
        xlabel(sprintf('%s percentile',valType_plot))
    end
    
end



