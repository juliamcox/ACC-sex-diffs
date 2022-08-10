function stats=latencyQuant_opto(cohort,behaviorTable,valType_plot,binNum,laserType,zscoreFlag,latencyType)


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
%% Extract animal means

groups = {'f';'m';'f_yfp';'m_yfp'};

for ng = 1:numel(groups)
    thisIDs = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIDs)
        for nb = 1:binNum
            eval(sprintf('X_%s(na,nb) = nanmean(behaviorTable.Latency(behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==0&behaviorTable.trial~=1));',groups{ng}))
            eval(sprintf('X_laser_%s(na,nb) = nanmean(behaviorTable.Latency(behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==laserType&behaviorTable.trial~=1));',groups{ng}))
            
            eval(sprintf('X_%s_count(na,nb) = nansum((behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==0&behaviorTable.trial~=1));',groups{ng}))
            eval(sprintf('X_laser_%s_count(na,nb) = nansum((behaviorTable.thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==laserType&behaviorTable.trial~=1));',groups{ng}))
        end
    end
end


%% Female laser v no laser by value
for nb = 1:binNum
        whichTest{1}(nb) = 1;
        [stats.pvals_posthoc(nb,1),~,stats.posthoc{nb,1}] = signrank(X_f(:,nb), X_laser_f(:,nb));
end

% Male laser v no laser by value
for nb = 1:binNum
        whichTest{2}(nb) = 1;
        [stats.pvals_posthoc(nb,2),~,stats.posthoc{nb,2}] = signrank(X_m(:,nb), X_laser_m(:,nb));
end

% female yfp
for nb = 1:binNum
        whichTest{3}(nb) = 1;
        [stats.pvals_posthoc(nb,3),~,stats.posthoc{nb,3}] = signrank(X_f_yfp(:,nb), X_laser_f_yfp(:,nb) );
end

% Male laser v no laser by value
for nb = 1:binNum
        whichTest{4}(nb) = 1;
        [stats.pvals_posthoc(nb,4),~,stats.posthoc{nb,4}] = signrank(X_m_yfp(:,nb), X_laser_m_yfp(:,nb));
end


%% Mixed-effects regession of  difference between laser and no laser

value = [];
LatencyDiff = [];
sex = [];
opsin = [];
animal = [];
for nb = 1:binNum
    LatencyDiff = cat(1,LatencyDiff,X_laser_f(:,nb)-X_f(:,nb),...
    X_laser_m(:,nb)-X_m(:,nb),...
    X_laser_f_yfp(:,nb)-X_f_yfp(:,nb),...
    X_laser_m_yfp(:,nb)-X_m_yfp(:,nb));

    value = cat(1,value,ones(size(X_laser_f(:,nb)-X_f(:,nb))).*nb,...
    ones(size(X_laser_m(:,nb)-X_m(:,nb))).*nb,...
    ones(size(X_laser_f_yfp(:,nb)-X_f_yfp(:,nb))).*nb,...
    ones(size(X_laser_m_yfp(:,nb)-X_m_yfp(:,nb))).*nb);
    sex = cat(1,sex,ones(size(X_laser_f(:,nb)-X_f(:,nb))),...
    zeros(size(X_laser_m(:,nb)-X_m(:,nb))),...
    ones(size(X_laser_f_yfp(:,nb)-X_f_yfp(:,nb))),...
    zeros(size(X_laser_m_yfp(:,nb)-X_m_yfp(:,nb))));
    opsin = cat(1,opsin,ones(size(X_laser_f(:,nb)-X_f(:,nb))),...
    ones(size(X_laser_m(:,nb)-X_m(:,nb))),...
    zeros(size(X_laser_f_yfp(:,nb)-X_f_yfp(:,nb))),...
    zeros(size(X_laser_m_yfp(:,nb)-X_m_yfp(:,nb))));
    
    animal = cat(1,animal,[1:size(X_f,1)]', [size(X_f,1)+1:size(X_f,1)+size(X_m,1)]',[size(X_f,1)+size(X_m)+1:size(X_f,1)+size(X_m,1)+size(X_f_yfp,1)]',[size(X_f,1)+size(X_m)+size(X_f_yfp)+1:size(X_f,1)+size(X_m,1)+size(X_f_yfp,1)+size(X_m_yfp,1)]');
    

end

T = table(LatencyDiff,'VariableNames',{'Latency'});
T.Sex = categorical(sex);
T.Opsin = categorical(opsin);
T.Value = categorical(value); 
T.Animal = categorical(animal);

f = 'Latency ~ Opsin*Value*Sex + (1+Value|Animal)';
glme_diff = fitlme(T,f,'DummyVarCoding','effects');
stats.glme_diff_coeff = dataset2cell((glme_diff.Coefficients));
stats.glme_diff=glme_diff;
stats.glme_diff_anova = dataset2cell(anova(glme_diff,'DFMethod','Satterthwaite'));

T = behaviorTable;
T = T(T.laser==0|T.laser==1,:);
T = T(~isnan(T.trialInit_thresh),:);
T.trialInit_thresh(T.trialInit_thresh==0) = 0.01;

T.latency = log(T.trialInit_thresh);
T.female = categorical(T.female);
T.qChosenDiff_quant = categorical(T.qChosenDiff_quant); 
T.opsin = categorical(T.opsin);
T.laser = categorical(T.laser);

f = 'latency ~ female*opsin*qChosenDiff_quant*laser + (1+qChosenDiff_quant|aID)';
mdl = fitlme(T,f);


%% Plot 
    
    
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;

   
  


figure();
% Plot female nphr outcome laser v no laser
subplot(2,2,1); hold on
plot(1:binNum,X_f,'Color',[.3 .3 .3 .7],'LineWidth',.25)
plot(1:binNum,X_laser_f,'Color',[plotParams.femaleC .7],'LineWidth',.25)
% 
% for nb = 1:binNum
%     xaxis = cat(2,repmat(nb,1,size(X_f,1))'-.2,repmat(nb,1,size(X_f,1))'+.2);
%     %plot(xaxis',cat(2,X_f(:,nb),X_laser_f(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_f,1))-.2;%+datasample(0:.01:.25,size(X_f,1));
%   % scatter(xaxis, X_f(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_f,1))+.2;%-datasample(0:.01:.25,size(X_f,1));
%   % scatter(xaxis, X_laser_f(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end

p=errorbar(nanmean(X_laser_f),nanstd(X_laser_f)./sqrt(size(X_laser_f,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(1:binNum,nanmean(X_f),nanstd(X_f)./sqrt(size(X_f,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
legend(p,{'laser';'no laser'})

subplot(2,2,2); hold on
plot(1:binNum,X_m,'Color',[.3 .3 .3 .7],'LineWidth',.25)
plot(1:binNum,X_laser_m,'Color',[plotParams.maleC .7],'LineWidth',.25)
% for nb = 1:binNum
%     xaxis = cat(2,repmat(nb,1,size(X_m,1))'-.2,repmat(nb,1,size(X_m,1))'+.2);
%     %plot(xaxis',cat(2,X_m(:,nb),X_laser_m(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_m,1))-.2;%+datasample(0:.01:.25,size(X_m,1));
%    %scatter(xaxis, X_m(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_m,1))+.2;%-datasample(0:.01:.25,size(X_m,1));
%   % scatter(xaxis, X_laser_m(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
p = errorbar(nanmean(X_m),nanstd(X_m)./sqrt(size(X_m,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(nanmean(X_laser_m),nanstd(X_laser_m)./sqrt(size(X_laser_m,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0);
legend(p,{'laser';'no laser'})

subplot(2,2,3); hold on
plot(1:binNum,X_f_yfp,'Color',[.3 .3 .3],'LineWidth',.25)
plot(1:binNum,X_laser_f_yfp,'Color',plotParams.femaleC,'LineWidth',.25)
% for nb = 1:binNum
%     xaxis = cat(2,repmat(nb,1,size(X_f_yfp,1))'-.2,repmat(nb,1,size(X_f_yfp,1))'+.2);
%    % plot(xaxis',cat(2,X_f_yfp(:,nb),X_laser_f_yfp(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
% end
errorbar(nanmean(X_f_yfp),nanstd(X_f_yfp)./sqrt(size(X_f_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(X_laser_f_yfp),nanstd(X_laser_f_yfp)./sqrt(size(X_laser_f_yfp,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0)
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_f_yfp,1))-.2;%+datasample(0:.01:.25,size(X_f_yfp,1));
%   % scatter(xaxis, X_f_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_f_yfp,1))+.2;%-datasample(0:.01:.25,size(X_f_yfp,1));
%  %  scatter(xaxis, X_laser_f_yfp(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end

subplot(2,2,4); hold on
plot(1:binNum,X_m_yfp,'Color',[.3 .3 .3],'LineWidth',.25)
plot(1:binNum,X_laser_m_yfp,'Color',plotParams.maleC,'LineWidth',.25)
% for nb = 1:binNum
%     xaxis = cat(2,repmat(nb,1,size(X_m_yfp,1))'-.2,repmat(nb,1,size(X_m_yfp,1))'+.2);
%   %  plot(xaxis',cat(2,X_m_yfp(:,nb),X_laser_m_yfp(:,nb))','-','Color',[.3 .3 .3 .8],'LineWidth',.25)
% end
errorbar(nanmean(X_m_yfp),nanstd(X_m_yfp)./sqrt(size(X_m_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(X_laser_m_yfp),nanstd(X_laser_m_yfp)./sqrt(size(X_laser_m_yfp,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0)
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_m_yfp,1))-.2;%+datasample(0:.01:.25,size(X_m_yfp,1));
%   % scatter(xaxis, X_m_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:binNum
%    xaxis = repmat(nb,1,size(X_m_yfp,1))+.2;%-datasample(0:.01:.25,size(X_m_yfp,1));
%  %  scatter(xaxis, X_laser_m_yfp(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end


for np = 1:4
    subplot(2,2,np)
    axis square
    ylabel('Trial initiation latency (sec)')
    set(gca,'YLim',[0,15],'XLim', [.5 max(unique(behaviorTable.thisValue_ptile))+.5])
    if laserType==3
        set(gca,'YLim',[0 15],'XLim', [.5 max(unique(behaviorTable.thisValue_ptile))+.5])
    end
    xlabel(sprintf('%s percentile',valType_plot))
end



