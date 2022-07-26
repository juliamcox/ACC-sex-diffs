function stats = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins)

%%% Figure 1

% Plot latency x value x sex as extracted in latencyXvalueXsex.m

% ext: file extension (LAS,UNI et)
% epochs: laser conditions to plot
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% latencyType: which latency to plot
% zscoreFlag: 0 or 1 to zscore resposnse times
% cohort: basename for cohort (e.g. 'ACC_DMS_nphr') for animal lists generated with generateAnimalList
% bins

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));




%% Plot trial initiation latency for males and females as a function of value (quantiles)
% Figure 1H
markersize = 10;
for nv = 1:numel(valType)
    for nl = 1:numel(latencyType)
        thisbin = eval(sprintf('bins.%s', valType{nv}));
        
        
        xaxis = 1:numel(thisbin);
        % female
        mu_f = eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv}));
        if isempty(mu_f)
            mu_f = nan(1,length(xaxis));
        end
        sem_f = nanstd(mu_f,[],1)./sqrt(size(mu_f,1));
        % male
        mu_m = eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv}));
        if isempty(mu_m)
            mu_m = nan(1,length(xaxis));
        end
        sem_m = nanstd(mu_m)./sqrt(size(mu_m,1));
        f=figure('Units','inches','Position',[5,5,2.5, 2.5]); hold on
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
            tempaxis = repmat(xaxis(np)+legendh(1).XOffset,1,size(mu_f,1))+datasample(-.04:.001:.04,size(mu_f,1));
            scatter(tempaxis,mu_f(:,np),10,'MarkerFaceColor',plotParams.femaleC,'MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        end
        
        for np = 1:numel(xaxis)
            tempaxis = repmat(xaxis(np)+legendh(2).XOffset,1,size(mu_m,1))+datasample(-.04:.001:.04,size(mu_m,1));
            scatter(tempaxis,mu_m(:,np),10,'MarkerFaceColor',plotParams.maleC,'MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        end
        
        set(gca,'FontSize',12,'XLim', [xaxis(1)-.5 xaxis(end-1)+.5])
        
        switch latencyType{nl}
            case 'trialStart'
                if zscoreFlag
                    ylabel('Trial initiation latency (zscore)')
                else
                    ylabel('Trial initiation latency (sec)')
                end
            case 'leverPress'
                if zscoreFlag
                    ylabel('Lever press latency (zscore)')
                else
                    ylabel('Lever press latency (sec)')
                end
            case 'withdraw'
                if zscoreFlag
                    ylabel('Nose poke withdrawal latency (zscore)')
                else
                    ylabel('Nose poke withdrawal latency (sec)')
                end
        end
        xlabel(sprintf('%s quantile',valType{nv}))
        l.Box = 'off';
        if zscoreFlag
            set(gca,'YLim',[-.3 .3])
        else
            set(gca,'YLim',[0 14.2])
        end
      
        if zscoreFlag
            set(gca,'YLim',[-.3 .3])
        else
            set(gca,'YLim',[0 14.2])
        end
        
        groupvec = [];
        vec = [];
        for na = 1:size(latency_f.trialInit_qChosenDiff,1)
            groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' na.*ones(size(thisbin(1:end-1)))'] );
            vec = cat(1,vec,mu_f(na,1:end-1)');
        end
        acounter = na+1;
        for na = 1:size(latency_m.trialInit_qChosenDiff,1)
            groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
            vec = cat(1,vec,mu_m(na,1:end-1)');
            acounter = acounter+1;
        end  
        % rank sum tests for male vs. female 
        for nb = 1:size(mu_f,2)-1
           [tempStats.postHoc_rankSum_p(nb),~,tempStats.postHoc_rankSum_stats{nb}] = ranksum(mu_f(:,nb),mu_m(:,nb)); 
        end
        eval(sprintf('stats.%s_%s_quant = tempStats;', valType{nv}, latencyType{nl}))
%         
% %         % mixed effects regression 
%         tbl = table(vec);
%         sex = categorical(groupvec(:,2));
%         
%         tbl.sex = sex;
%         tbl.value = categorical(groupvec(:,1),'ordinal',1);
%         subject = nominal(groupvec(:,3));
%         tbl.subject = subject;
%         
%         mdl = fitlme(tbl, 'vec~value+sex+value:sex+(1|subject)','DummyVarCoding','effects');
%         
%        % set up contrasts 
%         for nb = 1:max(groupvec(:,1))-1
%             h = zeros(1,numel(mdl.Coefficients.Name));
%             h(2) = 2;
%             h(strcmp(mdl.Coefficients.Name,sprintf('sex_0:value_%d',nb))) = 2;
%             [p(nb),fstat(nb),df1(nb),df2(nb)] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
%         end
%         h = zeros(1,numel(mdl.Coefficients.Name));
%         h(2) = 2;
%         h(end-nb+1:end) = -2;
%         [p(nb+1),fstat(nb+1),df1(nb+1),df2(nb+1)] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');        
%         eval(sprintf('stats.%s_%s_quant.coefTest.p = p;',valType{nv}, latencyType{nl}))
%         eval(sprintf('stats.%s_%s_quant.coefTest.fstat = fstat;',valType{nv}, latencyType{nl}))
%         eval(sprintf('stats.%s_%s_quant.coefTest.df1 = df1;',valType{nv}, latencyType{nl}))
%         eval(sprintf('stats.%s_%s_quant.coefTest.df2 = df2;',valType{nv}, latencyType{nl}))
%         eval(sprintf('stats.%s_%s_quant.mdl = mdl;',valType{nv}, latencyType{nl}))
%         
    end
end



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
    
    for na = 1:numel(thisplot)
        groupvec = cat(1,groupvec,[0,1,na+acounter]);
        vec = cat(1,vec,thisplot(na));
    end
    
    thisplot = eval(sprintf('latency_m.%s_prevNoReward;',latencyType{nl}));
    xaxis = repmat(b(2).XData(2)+b(2).XOffset,size(thisplot))-datasample(0:.001:.025,numel(thisplot))-.025;
    scatter(xaxis, thisplot, 12, 'MarkerEdgeColor','none', 'MarkerFaceColor', plotParams.maleC, 'MarkerFaceAlpha', .5)
    for na = 1:numel(thisplot)
        groupvec = cat(1,groupvec,[0,0,na+acounter]);
        vec = cat(1,vec,thisplot(na));
    end
    
    box off
    set(gca,'XTickLabel',{'Reward';'No reward'})
    xlabel('Previous trial')
    if zscoreFlag
        ylabel('Trial initiation latency (zscore)')
    else
        ylabel('Trial initiation latency (sec)')
        set(gca,'YLim',[0 14])
    end
    
    [stats.prevOutcome.p_nrew,~,stats.prevOutcome.stats_nrew] = ranksum(vec(groupvec(:,1)==0&groupvec(:,2)==0),vec(groupvec(:,1)==1&groupvec(:,2)==0));
    [stats.prevOutcome.p_rew,~,stats.prevOutcome.stats_rew] = ranksum(vec(groupvec(:,1)==0&groupvec(:,2)==1),vec(groupvec(:,1)==1&groupvec(:,2)==1));

     [stats.prevOutcome.p_f,~,stats.prevOutcome.stats_f] = signrank(vec(groupvec(:,1)==1&groupvec(:,2)==0),vec(groupvec(:,1)==1&groupvec(:,2)==1));
    [stats.prevOutcome.p_m,~,stats.prevOutcome.stats_m] = signrank(vec(groupvec(:,1)==0&groupvec(:,2)==0),vec(groupvec(:,1)==0&groupvec(:,2)==1));

    tbl = table(vec);
    sex = categorical(groupvec(:,1));
    
    tbl.sex = sex;
    tbl.outcome = categorical(groupvec(:,2),'ordinal',1);
    subject = nominal(groupvec(:,3));
    tbl.subject = subject;
    
    mdl = fitlme(tbl, 'vec~outcome+sex+outcome:sex+(1|subject)','DummyVarCoding','effects');
    stats.prevOutcome.mdl = mdl;
    stats.prevOutcome.mdl_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite'));
    stats.prevOutcome.coeff = dataset2cell(mdl.Coefficients); 
    

    h_m = [1 1 1 1];
    h_f = [1 -1 1 -1];
    h = h_m-h_f;
    [ stats.prevOutcome.constrasts.p_nrew, stats.prevOutcome.constrasts.fstat_nrew, stats.prevOutcome.constrasts.df1_nrew, stats.prevOutcome.constrasts.df2_nrew] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
   
    h_m = [1 1 -1 -1];
    h_f = [1 -1 -1 1];
    h = h_m-h_f;
    [ stats.prevOutcome.constrasts.p_rew, stats.prevOutcome.constrasts.fstat_rew, stats.prevOutcome.constrasts.df1_rew, stats.prevOutcome.constrasts.df2_rew] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');

    h_r = [1 -1 1 -1];
    h_nr = [1 -1 -1 1];
    h = h_r-h_nr;
    [ stats.prevOutcome.constrasts.p_f, stats.prevOutcome.constrasts.fstat_f, stats.prevOutcome.constrasts.df1_f, stats.prevOutcome.constrasts.df2_f] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');
    
    h_r = [1 1 1 1];
    h_nr = [1 1 -1 -1];
    h = h_r-h_nr;    [ stats.prevOutcome.constrasts.p_m, stats.prevOutcome.constrasts.fstat_m, stats.prevOutcome.constrasts.df1_m, stats.prevOutcome.constrasts.df2_m] = coefTest(mdl,h,zeros(size(h,1)),'DFMethod','Satterthwaite');

end