function stats = latencyXoutcomeXsex_plot(latency_f,latency_m,latencyType,zscoreFlag,behaviorTable,aids_f,aids_m)

% latency_f and latency_m: structures with animal averages of latency x value and previous outcome x value
% latencyType: trialInit or lever press 
% zscoreFlag: zscored latencies?

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

%% Plot latency as a function of previous outcome
% Figure 1C
for nl = 1:numel(latencyType)
    groupvec = [];
    vec = [];
    
    f=figure('Units','inches','Position',[5,5,1.8 2.5]);
    
    mu = [eval(sprintf('mean(latency_f.%s_prevReward)',latencyType{nl})) eval(sprintf('mean(latency_f.%s_prevNoReward)',latencyType{nl}));eval(sprintf('mean(latency_m.%s_prevReward)',latencyType{nl})) eval(sprintf('mean(latency_m.%s_prevNoReward)',latencyType{nl}))];
    sem = [eval(sprintf('std(latency_f.%s_prevReward)./sqrt(length(latency_f.%s_prevReward))',latencyType{nl},latencyType{nl})) eval(sprintf('std(latency_f.%s_prevNoReward)./sqrt(length(latency_f.%s_prevReward))',latencyType{nl},latencyType{nl}));eval(sprintf('std(latency_m.%s_prevReward)./sqrt(length(latency_m.%s_prevReward))',latencyType{nl},latencyType{nl})) eval(sprintf('std(latency_m.%s_prevNoReward)./sqrt(length(latency_m.%s_prevReward))',latencyType{nl},latencyType{nl}))];
    b=bar(mu');
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
    thisplot = eval(sprintf('latency_f.%s_prevReward;',latencyType{nl}));
    xaxis = repmat(b(1).XData(1)+b(1).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
    scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)
    % Set up data for mixed effects regression
    for na = 1:numel(thisplot)
        groupvec = cat(1,groupvec,[1,1,na]); %sex, outcome, animalID)
        vec = cat(1,vec,thisplot(na));
    end
    
    thisplot = eval(sprintf('latency_f.%s_prevNoReward;',latencyType{nl}));
    xaxis = repmat(b(1).XData(2)+b(1).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
    scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)
    for na = 1:numel(thisplot)
        groupvec = cat(1,groupvec,[1,0,na]);
        vec = cat(1,vec,thisplot(na));
    end
    acounter = na;
    thisplot = eval(sprintf('latency_m.%s_prevReward;',latencyType{nl}));
    xaxis = repmat(b(2).XData(1)+b(2).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
    scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
    
  
    thisplot = eval(sprintf('latency_m.%s_prevNoReward;',latencyType{nl}));
    xaxis = repmat(b(2).XData(2)+b(2).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
    scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
    
    box off
    set(gca,'XTickLabel',{'Reward';'No reward'})
    xlabel('Previous trial')
    if zscoreFlag
        ylabel('Trial initiation latency (zscore)')
    else
        ylabel('Trial initiation latency (sec)')
        set(gca,'YLim',[0 14])
    end
   
end


for nl = 1:numel(latencyType)
    % select animals
    idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
    idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
    X = behaviorTable(cell2mat(idx),:);
    % select non-laser sessions
    X = X(X.laserSession==0,:);
    X.trial = nanzscore(X.trial);
    X.female = categorical(X.female);
    X.previousReward = categorical(X.previousReward);
    f = sprintf('%s ~ female*previousReward + (1+previousReward|aID)',latencyType{nl});
    mdl = fitlme(X,f,'DummyVarCoding','effects');
    temp.mdl = mdl;
    temp.mdl_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite'));
    temp.mdl_coeff = dataset2cell(mdl.Coefficients);
    
    h_m_rew = [1 -1 1 -1];
    h_m_nrew = [1 1 1 1];
    h_f_rew = [1 -1 -1 1];
    h_f_nrew = [1 1 -1 -1];
    h = h_m_rew-h_f_rew;
    [temp.contrasts.rew_p_p,temp.contrasts.rew_p_f,temp.contrasts.rew_p_df1,temp.contrasts.rew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
    h = h_m_nrew-h_f_nrew;
    [temp.contrasts.nrew_p_p,temp.contrasts.nrew_p_f,temp.contrasts.nrew_p_df1,temp.contrasts.nrew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
    h = h_f_rew-h_f_nrew;
    [temp.contrasts.f_p_p,temp.contrasts.f_p_f,temp.contrasts.f_p_df1,temp.contrasts.f_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
    h = h_m_rew-h_m_nrew;
    [temp.contrasts.m_p_p,temp.contrasts.m_p_f,temp.contrasts.m_p_df1,temp.contrasts.m_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
    eval(sprintf('stats.%s = temp;', latencyType{nl}));
end

%% 
% 
% %% Plot latency as a function of previous outcome
% % Figure 1C
% for nl = 1:numel(latencyType)
%     groupvec = [];
%     vec = [];
%     
%     f=figure('Units','inches','Position',[5,5,1.8 2.5]);
%     
%     mu = [eval(sprintf('mean(latency_f.%s_prevReward)',latencyType{nl})) eval(sprintf('mean(latency_f.%s_prevNoReward)',latencyType{nl}));eval(sprintf('mean(latency_m.%s_prevReward)',latencyType{nl})) eval(sprintf('mean(latency_m.%s_prevNoReward)',latencyType{nl}))];
%     sem = [eval(sprintf('std(latency_f.%s_prevReward)./sqrt(length(latency_f.%s_prevReward))',latencyType{nl},latencyType{nl})) eval(sprintf('std(latency_f.%s_prevNoReward)./sqrt(length(latency_f.%s_prevReward))',latencyType{nl},latencyType{nl}));eval(sprintf('std(latency_m.%s_prevReward)./sqrt(length(latency_m.%s_prevReward))',latencyType{nl},latencyType{nl})) eval(sprintf('std(latency_m.%s_prevNoReward)./sqrt(length(latency_m.%s_prevReward))',latencyType{nl},latencyType{nl}))];
%     b=bar(mu');
%     b(1).FaceColor = 'none';
%     b(1).EdgeColor = plotParams.femaleC;
%     b(2).FaceColor = 'none';
%     b(2).EdgeColor = plotParams.maleC;
%     b(1).LineWidth = 1.5;
%     b(2).LineWidth = 1.5;
%     pause(.001)
%     hold on
%     errorbar(b(1).XData+b(1).XOffset, b(1).YData,sem(1,:),'LineStyle','none','Color',plotParams.femaleC,'LineWidth',1.5,'CapSize',0);
%     errorbar(b(2).XData+b(2).XOffset, b(2).YData,sem(2,:),'LineStyle','none','Color',plotParams.maleC,'LineWidth',1.5,'CapSize',0);
%     %Plot individual data points
%     thisplot = eval(sprintf('latency_f.%s_prevReward;',latencyType{nl}));
%     xaxis = repmat(b(1).XData(1)+b(1).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
%     scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)
%     % Set up data for mixed effects regression
%     for na = 1:numel(thisplot)
%         groupvec = cat(1,groupvec,[1,1,na]); %sex, outcome, animalID)
%         vec = cat(1,vec,thisplot(na));
%     end
%     
%     thisplot = eval(sprintf('latency_f.%s_prevNoReward;',latencyType{nl}));
%     xaxis = repmat(b(1).XData(2)+b(1).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
%     scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.femaleC, 'MarkerFaceAlpha', .5)
%     for na = 1:numel(thisplot)
%         groupvec = cat(1,groupvec,[1,0,na]);
%         vec = cat(1,vec,thisplot(na));
%     end
%     acounter = na;
%     thisplot = eval(sprintf('latency_m.%s_prevReward;',latencyType{nl}));
%     xaxis = repmat(b(2).XData(1)+b(2).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
%     scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
%     
%     for na = 1:numel(thisplot)
%         groupvec = cat(1,groupvec,[0,1,na+acounter]);
%         vec = cat(1,vec,thisplot(na));
%     end
%     
%     thisplot = eval(sprintf('latency_m.%s_prevNoReward;',latencyType{nl}));
%     xaxis = repmat(b(2).XData(2)+b(2).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
%     scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
%     for na = 1:numel(thisplot)
%         groupvec = cat(1,groupvec,[0,0,na+acounter]);
%         vec = cat(1,vec,thisplot(na));
%     end
%     
%     box off
%     set(gca,'XTickLabel',{'Reward';'No reward'})
%     xlabel('Previous trial')
%     if zscoreFlag
%         ylabel('Trial initiation latency (zscore)')
%     else
%         ylabel('Trial initiation latency (sec)')
%         set(gca,'YLim',[0 14])
%     end
%     
%     [temp.prevOutcome.p_nrew,~,temp.prevOutcome.temp_nrew] = ranksum(vec(groupvec(:,1)==0&groupvec(:,2)==0),vec(groupvec(:,1)==1&groupvec(:,2)==0));
%     [temp.prevOutcome.p_rew,~,temp.prevOutcome.temp_rew] = ranksum(vec(groupvec(:,1)==0&groupvec(:,2)==1),vec(groupvec(:,1)==1&groupvec(:,2)==1));
% 
%      [temp.prevOutcome.p_f,~,temp.prevOutcome.temp_f] = signrank(vec(groupvec(:,1)==1&groupvec(:,2)==0),vec(groupvec(:,1)==1&groupvec(:,2)==1));
%     [temp.prevOutcome.p_m,~,temp.prevOutcome.temp_m] = signrank(vec(groupvec(:,1)==0&groupvec(:,2)==0),vec(groupvec(:,1)==0&groupvec(:,2)==1));
% 
%     tbl = table(vec);
%     sex = categorical(groupvec(:,1));
%     
%     tbl.sex = sex;
%     tbl.outcome = categorical(groupvec(:,2),'ordinal',1);
%     subject = nominal(groupvec(:,3));
%     tbl.subject = subject;
%     
%     mdl = fitlme(tbl, 'vec~outcome+sex+outcome:sex+(1+outcome|subject)','DummyVarCoding','effects');
%     temp.prevOutcome.mdl = mdl;
%     temp.prevOutcome.mdl_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite'));
%     temp.prevOutcome.coeff = dataset2cell(mdl.Coefficients); 
%     
% 
%     h_m = [1 1 1 1];
%     h_f = [1 -1 1 -1];
%     h = h_m-h_f;
%     [ temp.prevOutcome.constrasts.p_nrew, temp.prevOutcome.constrasts.fstat_nrew, temp.prevOutcome.constrasts.df1_nrew, temp.prevOutcome.constrasts.df2_nrew] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
%    
%     h_m = [1 1 -1 -1];
%     h_f = [1 -1 -1 1];
%     h = h_m-h_f;
%     [ temp.prevOutcome.constrasts.p_rew, temp.prevOutcome.constrasts.fstat_rew, temp.prevOutcome.constrasts.df1_rew, temp.prevOutcome.constrasts.df2_rew] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
% 
%     h_r = [1 -1 1 -1];
%     h_nr = [1 -1 -1 1];
%     h = h_r-h_nr;
%     [ temp.prevOutcome.constrasts.p_f, temp.prevOutcome.constrasts.fstat_f, temp.prevOutcome.constrasts.df1_f, temp.prevOutcome.constrasts.df2_f] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
%     
%     h_r = [1 1 1 1];
%     h_nr = [1 1 -1 -1];
%     h = h_r-h_nr;    [ temp.prevOutcome.constrasts.p_m, temp.prevOutcome.constrasts.fstat_m, temp.prevOutcome.constrasts.df1_m, temp.prevOutcome.constrasts.df2_m] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
% 
% end