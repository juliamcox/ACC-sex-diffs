function stats = choiceXvalueXsex_plot(choice_f, choice_m, valType, bins,behaviorTable, aids_f,aids_m)

%%% Figure 1B,G

% Plot stay prob and probability right choice x value x sex as extracted in ctrl_valueXsex.m

% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% bins 

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));



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

%% Mixed-effects regression
behaviorTable = behaviorTable(behaviorTable.laserSession == 0,:);
for nv = 1:numel(valType)
%% Mixed-effects model with weight, value, trial number and sex 
% select animals
idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
X = behaviorTable(cell2mat(idx),:);
% select non-laser sessions
X = X(X.laserSession==0,:);
X = X(~isnan(X.trialInit_thresh),:);
X.trial = nanzscore(X.trial); 
X.female = categorical(X.female);
eval(sprintf('X.%s_quant_choice = categorical(X.%s_quant_choice);',valType{nv},valType{nv}));

f = sprintf('choice ~ %s_quant_choice*female*trial + (1+%s_quant_choice*trial|aID)',valType{nv},valType{nv});

mdl = fitglme(X,f,'DummyVarCoding','effects','Distribution','binomial');
eval(sprintf('stats.mdl_%s_quant = mdl;', valType{nv})); 
eval(sprintf('stats.mdl_anova_%s_quant = dataset2cell(anova(mdl));',valType{nv})); 
eval(sprintf('stats.mdl_coeff_%s_quant = dataset2cell(mdl.Coefficients);',valType{nv})); 


end


