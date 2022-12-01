function stats=choice_opto(behaviorTable,valType_plot,binNum,laserType,cohort)

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
behaviorTable(isnan(behaviorTable.trialInit_thresh),:) = [];

thisValue_ptile = eval(sprintf('behaviorTable.%s_quant_choice',valType_plot));% which type of value
thisValue = eval(sprintf('behaviorTable.%s',valType_plot)); % which type of value


%% Extract choice data for each animal 

ids_m = generateAnimalList(sprintf('%s_male',cohort));
ids_f = generateAnimalList(sprintf('%s_female',cohort));
ids_m_yfp = generateAnimalList(sprintf('%s_yfp_male',cohort));
ids_f_yfp = generateAnimalList(sprintf('%s_yfp_female',cohort));

groups = {'f';'m';'f_yfp';'m_yfp'};

for ng = 1:numel(groups)
    thisIDs = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIDs)
        for nb = 1:binNum
            eval(sprintf('X_%s(na,nb) = nanmean(behaviorTable.choice(thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==0&behaviorTable.laserSession==1)==0);',groups{ng}))
            eval(sprintf('X_laser_%s(na,nb) = nanmean(behaviorTable.choice(thisValue_ptile==nb&strcmp(behaviorTable.aID,thisIDs{na})&behaviorTable.laser==laserType&behaviorTable.laserSession==1)==0);',groups{ng}))
        end
    end
end



%%
%% Mixed-effects regession of  difference between laser and no laser

value = [];
ChoiceDiff = [];
sex = [];
opsin = [];
animal = [];
for nb = 1:binNum
    ChoiceDiff = cat(1,ChoiceDiff,X_laser_f(:,nb)-X_f(:,nb),...
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

T = table(ChoiceDiff,'VariableNames',{'Choice'});
T.Sex = categorical(sex);
T.Opsin = categorical(opsin);
T.Value = categorical(value); 
T.Animal = categorical(animal);

f = 'Choice ~ Opsin*Value*Sex + (1+Value|Animal)';
glme_diff = fitlme(T,f,'DummyVarCoding','effects');
stats.glme_diff_coeff = dataset2cell((glme_diff.Coefficients));
stats.glme_diff=glme_diff;
stats.glme_diff_anova = dataset2cell(anova(glme_diff,'DFMethod','Satterthwaite'));

%% Female laser v no laser by value
for nb = 1:binNum
        whichTest{1}(nb) = 1;
        [stats.pvals_posthoc_choice(nb,1),~,stats.posthoc_choice{nb,1}] = signrank(X_f(:,nb), X_laser_f(:,nb),'Method','approximate');
end

% Male laser v no laser by value
for nb = 1:binNum

        whichTest{2}(nb) = 1;
        [stats.pvals_posthoc_choice(nb,2),~,stats.posthoc_choice{nb,2}] = signrank(X_m(:,nb), X_laser_m(:,nb),'Method','approximate');
end

% female yfp
for nb = 1:binNum
        whichTest{3}(nb) = 1;
        [stats.pvals_posthoc_choice(nb,3),~,stats.posthoc_choice{nb,3}] = signrank(X_f_yfp(:,nb), X_laser_f_yfp(:,nb),'Method','approximate' );
end

% Male laser v no laser by value
for nb = 1:binNum
        whichTest{4}(nb) = 1;
        [stats.pvals_posthoc_choice(nb,4),~,stats.posthoc_choice{nb,4}] = signrank(X_m_yfp(:,nb), X_laser_m_yfp(:,nb),'Method','approximate');
end



%% plot choice
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;
clear mu_f mu_f_laser mu_f_yfp mu_m mu_m_laser mu_m_yfp mu_m_yfp_laser mu_f_yfp_laser
ids_f = unique(behaviorTable.aID(behaviorTable.female==1&behaviorTable.opsin==1));
ids_f_yfp = unique(behaviorTable.aID(behaviorTable.female==1&behaviorTable.opsin==0));
ids_m = unique(behaviorTable.aID(behaviorTable.female==0&behaviorTable.opsin==1));
ids_m_yfp = unique(behaviorTable.aID(behaviorTable.female==0&behaviorTable.opsin==0));

trialNums = arrayfun(@(x) numel(behaviorTable.choice(strcmp(behaviorTable.aID,thisIDs{na})  & behaviorTable.laser==laserType&behaviorTable.laserSession==1)), ids_f);
ids_f = ids_f(trialNums>0);
trialNums = arrayfun(@(x) numel(behaviorTable.choice(strcmp(behaviorTable.aID,thisIDs{na})  & behaviorTable.laser==laserType&behaviorTable.laserSession==1)), ids_f_yfp);
ids_f_yfp = ids_f_yfp(trialNums>0);

trialNums = arrayfun(@(x) numel(behaviorTable.choice(strcmp(behaviorTable.aID,thisIDs{na})  & behaviorTable.laser==laserType&behaviorTable.laserSession==1)), ids_m);
ids_m = ids_m(trialNums>0);
trialNums = arrayfun(@(x) numel(behaviorTable.choice(strcmp(behaviorTable.aID,thisIDs{na})  & behaviorTable.laser==laserType&behaviorTable.laserSession==1)), ids_m_yfp);
ids_m_yfp = ids_m_yfp(trialNums>0);



for nb = 1:numel(unique(thisValue_ptile))
    for na = 1:numel(ids_f)
        mu_f(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_f{na})  & behaviorTable.laser==0 & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
        mu_f_laser(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_f{na})  & behaviorTable.laser==laserType & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
    end
    for na = 1:numel(ids_f_yfp)
        mu_f_yfp(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_f_yfp{na})  & behaviorTable.laser==0 & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
        mu_f_yfp_laser(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_f_yfp{na})  & behaviorTable.laser==laserType & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
    end
    for na = 1:numel(ids_m)
        mu_m(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_m{na})  & behaviorTable.laser==0 & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
        mu_m_laser(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_m{na})  & behaviorTable.laser==laserType & thisValue_ptile==nb & behaviorTable.laserSession == 1  )==0));
    end
    for na = 1:numel(ids_m_yfp)
        mu_m_yfp(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_m_yfp{na})  & behaviorTable.laser==0 & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
        mu_m_yfp_laser(na,nb) = nanmean((behaviorTable.choice(strcmp(behaviorTable.aID,ids_m_yfp{na}) & behaviorTable.laser==laserType & thisValue_ptile==nb & behaviorTable.laserSession == 1 )==0));
    end
end



figure();
% Plot female nphr outcome laser v no laser
subplot(3,2,1); hold on
% 
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_f,1));%+datasample(0:.01:.25,size(mu_f,1));
%  %  scatter(xaxis, mu_f(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_f,1));%-datasample(0:.01:.25,size(mu_f,1));
%  %  scatter(xaxis, mu_f_laser(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
p=errorbar(nanmean(mu_f_laser),nanstd(mu_f_laser)./sqrt(size(mu_f_laser,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(1:numel(unique(thisValue_ptile)),nanmean(mu_f),nanstd(mu_f)./sqrt(size(mu_f,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);


subplot(3,2,3); hold on
plot(1:binNum,mu_f,'Color',[.3 .3 .3 .7],'LineWidth',1)
plot(1:binNum, mu_f_laser,'Color',[plotParams.femaleC .7],'LineWidth',1)


subplot(3,2,2); hold on

% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_m,1));%+datasample(0:.01:.25,size(mu_m,1));
%  %  scatter(xaxis, mu_m(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_m,1));%-datasample(0:.01:.25,size(mu_m,1));
%   % scatter(xaxis, mu_m_laser(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
p = errorbar(nanmean(mu_m),nanstd(mu_m)./sqrt(size(mu_m,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(nanmean(mu_m_laser),nanstd(mu_m_laser)./sqrt(size(mu_m_laser,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0);


%legend(p,{'laser';'no laser'})
subplot(3,2,4); hold on

plot(1:binNum,mu_m','Color',[.3 .3 .3 .7],'LineWidth',1)
plot(1:binNum, mu_m_laser','Color',[plotParams.maleC .7],'LineWidth',1)




subplot(3,2,5); hold on

plot(1:binNum,mu_f_yfp,'Color',[.3 .3 .3 .7],'LineWidth',1)
plot(1:binNum, mu_f_yfp_laser,'Color',[plotParams.femaleC .7],'LineWidth',1)
errorbar(nanmean(mu_f_yfp),nanstd(mu_f_yfp)./sqrt(size(mu_f_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(mu_f_yfp_laser),nanstd(mu_f_yfp_laser)./sqrt(size(mu_f_yfp_laser,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0)
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_f_yfp,1));%+datasample(0:.01:.25,size(mu_f_yfp,1));
%   % scatter(xaxis, mu_f_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_f_yfp,1));%-datasample(0:.01:.25,size(mu_f_yfp,1));
%   % scatter(xaxis, mu_f_yfp_laser(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end

subplot(3,2,6); hold on

plot(1:binNum,mu_m_yfp,'Color',[.3 .3 .3 .7],'LineWidth',1)
plot(1:binNum, mu_m_yfp_laser,'Color',[plotParams.maleC .7],'LineWidth',1)
errorbar(nanmean(mu_m_yfp),nanstd(mu_m_yfp)./sqrt(size(mu_m_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(mu_m_yfp_laser),nanstd(mu_m_yfp_laser)./sqrt(size(mu_m_yfp_laser,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0)
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_m_yfp,1));%+datasample(0:.01:.25,size(mu_m_yfp,1));
%   % scatter(xaxis, mu_m_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end
% for nb = 1:numel(unique(thisValue_ptile))
%    xaxis = repmat(nb,1,size(mu_m_yfp,1));%-datasample(0:.01:.25,size(mu_m_yfp,1));
%  %  scatter(xaxis, mu_m_yfp_laser(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
% end


for np = 1:6
    subplot(3,2,np)
    ylabel('P(choice=R)')
    set(gca,'YLim',[0,1],'XLim', [.5 max(unique(thisValue_ptile))+.5])
    if laserType==3
        set(gca,'YLim',[0 1],'XLim', [.5 max(unique(thisValue_ptile))+.5])
    end
    xlabel('thisValue percentile')
end