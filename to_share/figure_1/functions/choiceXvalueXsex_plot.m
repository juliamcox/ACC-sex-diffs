function choiceXvalueXsex_plot(choice_f, choice_m, valType, bins)

%%% Figure 1B,G

% Plot stay prob and probability right choice x value x sex as extracted in ctrl_valueXsex.m

% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% bins 

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); 



%% Plot probability right choice for males and females as a function of value quantiles
% Figure 1G
markersize = nan;
for nv = 1:numel(valType)
    
    f=figure('Units','inches','Position',[5,5,1.8 2.5]);
    thisbin = eval(sprintf('bins.%s', valType{nv}));
    
    % Plot female data
    hold on
    xaxis = 1:numel(thisbin);
    
    mu_f = eval(sprintf('choice_f.choice_%s_quant;',valType{nv}));
    if isempty(mu_f)
        mu_f = nan(1,length(xaxis));
    end
    sem_f = nanstd(mu_f,[],1)./sqrt(size(mu_f,1));
    
    mu_m = eval(sprintf('choice_m.choice_%s_quant;',valType{nv}));
    if isempty(mu_m)
        mu_m = nan(1,length(xaxis));
    end
    sem_m = nanstd(mu_m)./sqrt(size(mu_m,1));
    
    
    legendh = bar(xaxis, [nanmean(mu_f,1); nanmean(mu_m,1)]');
    legendh(1).EdgeColor = 'none';
    legendh(2).EdgeColor = 'none';
    legendh(1).FaceColor = 'none';
    legendh(2).FaceColor = 'none';
    legendh(1).LineWidth = 1.5;
    legendh(2).LineWidth = 1.5;
    
    pause(.001)
    
    errorbar(xaxis,nanmean(mu_f,1), sem_f, 'LineStyle', '-', 'CapSize', 0, 'Color', plotParams.femaleC, 'LineWidth', 1);
    errorbar(xaxis,nanmean(mu_m,1), sem_m, 'LineStyle', '-', 'CapSize', 0, 'Color', plotParams.maleC, 'LineWidth', 1);
    
    for np = 1:numel(xaxis)
        tempaxis = repmat(xaxis(np)+legendh(1).XOffset,1,size(mu_f,1));%+datasample(-.04:.001:.04,size(mu_f,1));
        scatter(tempaxis,mu_f(:,np),10,'MarkerFaceColor',plotParams.femaleC,'MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
    end
    
    for np = 1:numel(xaxis)
        tempaxis = repmat(xaxis(np)+legendh(2).XOffset,1,size(mu_m,1));%+datasample(-.04:.001:.04,size(mu_m,1));
        scatter(tempaxis,mu_m(:,np),10,'MarkerFaceColor',plotParams.maleC,'MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
    end
    
    
    set(gca,'FontSize',12,'XLim', [xaxis(1)-1 xaxis(end)])
    ylabel('P(choice==right)')
    xlabel(sprintf('%s quantile bin',valType{nv}));
 
end




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



end


