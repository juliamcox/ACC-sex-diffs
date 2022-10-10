function stats = choiceXoutcomeXsex_plot(choice_f, choice_m,behaviorTable,aids_f,aids_m)

%%% Figure 1b

% Plot stay probability by previous outcome

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));


%% Plot stay probability by previous outcome

mu  = [mean(choice_f.stayProb_prevReward) mean(choice_f.stayProb_prevNoReward); mean(choice_m.stayProb_prevReward) mean(choice_m.stayProb_prevNoReward)];
sem = [std(choice_f.stayProb_prevReward)./sqrt(length(choice_f.stayProb_prevReward)) std(choice_f.stayProb_prevNoReward)./sqrt(length(choice_f.stayProb_prevReward)); std(choice_m.stayProb_prevReward)./sqrt(length(choice_m.stayProb_prevReward)) std(choice_m.stayProb_prevNoReward)./sqrt(length(choice_m.stayProb_prevReward))];


f=figure('Units','inches','Position',[5,5,1.8 2.5]);

b = bar(mu');
b(1).FaceColor = 'none';
b(1).EdgeColor = plotParams.femaleC;
b(2).FaceColor = 'none';
b(2).EdgeColor = plotParams.maleC;
b(1).LineWidth = 1.5;
b(2).LineWidth = 1.5; 
pause(.001)
hold on
errorbar(b(1).XData+b(1).XOffset, b(1).YData,sem(1,:),'LineStyle','none','Color',plotParams.femaleC,'LineWidth',1.5,'CapSize',0);
errorbar(b(2).XData+b(2).XOffset, b(2).YData,sem(2,:),'LineStyle','none','Color',plotParams.maleC,'LineWidth',1.5,'CapSize',0);
%Plot individual data points
xaxis = repmat(b(1).XData(1)+b(1).XOffset,size(choice_f.stayProb_prevReward))-datasample(0:.001:.025,numel(choice_f.stayProb_prevReward))-.025;
scatter(xaxis, choice_f.stayProb_prevReward, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)

xaxis = repmat(b(1).XData(2)+b(1).XOffset,size(choice_f.stayProb_prevNoReward))-datasample(0:.001:.025,numel(choice_f.stayProb_prevNoReward))-.025;
scatter(xaxis, choice_f.stayProb_prevNoReward, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)
xaxis = repmat(b(2).XData(1)+b(2).XOffset,size(choice_m.stayProb_prevReward))-datasample(0:.001:.025,numel(choice_m.stayProb_prevReward))-.025;
scatter(xaxis, choice_m.stayProb_prevReward, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
xaxis = repmat(b(2).XData(2)+b(2).XOffset,size(choice_m.stayProb_prevNoReward))-datasample(0:.001:.025,numel(choice_m.stayProb_prevNoReward))-.025;
scatter(xaxis, choice_m.stayProb_prevNoReward, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)

set(gca,'YLim',[0 1], 'XTickLabel', {'Reward';'No reward'})
xlabel('Previous outcome')
ylabel('Stay probability') 
box off

%% Mixed-effects regression
behaviorTable = behaviorTable(behaviorTable.laserSession == 0,:);
% select animals
idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
X = behaviorTable(cell2mat(idx),:);
% select non-laser sessions
X = X(X.laserSession==0,:);
X = X(~isnan(X.trialInit_thresh),:);
X.trial = nanzscore(X.trial); 
X.female = categorical(X.female);
X.previousReward = categorical(X.previousReward); 

f = 'stay ~ previousReward*female*trial + (1+previousReward*trial|aID)';

mdl = fitglme(X,f,'DummyVarCoding','effects','Distribution','Binomial');
stats.mdl = mdl;
stats.mdl_anova = dataset2cell(anova(mdl)); 
stats.mdl_coeff = dataset2cell(mdl.Coefficients);

end


