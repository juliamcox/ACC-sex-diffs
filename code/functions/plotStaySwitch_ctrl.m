function stats = plotStaySwitch_ctrl(latency_f,latency_m,latencyType)

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
nf = numel(latency_f.trialInit_switch);
nm = numel(latency_m.trialInit_switch); 

for nl = 1:numel(latencyType)
    mu = [mean(eval(sprintf('latency_f.%s_stay',latencyType{nl}))) mean(eval(sprintf('latency_f.%s_switch',latencyType{nl})));...
        mean(eval(sprintf('latency_m.%s_stay',latencyType{nl}))) mean(eval(sprintf('latency_m.%s_switch',latencyType{nl})))];
    
    sem = [std(eval(sprintf('latency_f.%s_stay',latencyType{nl})))./sqrt(nf) std(eval(sprintf('latency_f.%s_switch',latencyType{nl})))./sqrt(nf);...
        std(eval(sprintf('latency_m.%s_stay',latencyType{nl}))./sqrt(nm)) std(eval(sprintf('latency_m.%s_switch',latencyType{nl})))./sqrt(nm)];
    
    f=figure('Units','inches','Position',[5,5,7 2.5]); 
    subplot(1,3,1);hold on
    b=bar(mu');
    b(1).EdgeColor = plotParams.femaleC;
    b(1).FaceColor = 'none';
    b(2).EdgeColor = plotParams.maleC;
    b(2).FaceColor = 'none';
    b(2).LineWidth = 1.5;
    b(1).LineWidth = 1.5;
    pause(.001);
    errorbar(b(1).XData+b(1).XOffset,mu(1,:),sem(1,:), 'Color',plotParams.femaleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    errorbar(b(2).XData+b(2).XOffset,mu(2,:),sem(2,:), 'Color',plotParams.maleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    scatter(repmat(b(1).XData(1)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_stay',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(1).XData(2)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_switch',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(1)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_stay',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(2)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_switch',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)

    set(gca,'XTick',b(1).XData,'XTickLabel',{'stay';'switch'},'YLim',[0 15])
    ylabel(latencyType{nl})
    
    mu = [mean(eval(sprintf('latency_f.%s_stayRew',latencyType{nl}))) mean(eval(sprintf('latency_f.%s_switchRew',latencyType{nl})));...
        mean(eval(sprintf('latency_m.%s_stayRew',latencyType{nl}))) mean(eval(sprintf('latency_m.%s_switchRew',latencyType{nl})))];
    
    sem = [std(eval(sprintf('latency_f.%s_stayRew',latencyType{nl})))./sqrt(nf) std(eval(sprintf('latency_f.%s_switchRew',latencyType{nl})))./sqrt(nf);...
        std(eval(sprintf('latency_m.%s_stayRew',latencyType{nl}))./sqrt(nm)) std(eval(sprintf('latency_m.%s_switchRew',latencyType{nl})))./sqrt(nm)];
    
    subplot(1,3,2); hold on
    b=bar(mu');
    b(1).EdgeColor = plotParams.femaleC;
    b(1).FaceColor = 'none';
    b(2).EdgeColor = plotParams.maleC;
    b(2).FaceColor = 'none';
    b(2).LineWidth = 1.5;
    b(1).LineWidth = 1.5;
    pause(.001);
    errorbar(b(1).XData+b(1).XOffset,mu(1,:),sem(1,:), 'Color',plotParams.femaleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    errorbar(b(2).XData+b(2).XOffset,mu(2,:),sem(2,:), 'Color',plotParams.maleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    scatter(repmat(b(1).XData(1)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_stayRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(1).XData(2)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_switchRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(1)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_stayRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(2)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_switchRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)

    set(gca,'XTick',b(1).XData,'XTickLabel',{'stay';'switch'},'YLim',[0 15])
    ylabel(latencyType{nl})
    title('Previously rewarded trials')
    
    mu = [mean(eval(sprintf('latency_f.%s_stayNRew',latencyType{nl}))) mean(eval(sprintf('latency_f.%s_switchNRew',latencyType{nl})));...
        mean(eval(sprintf('latency_m.%s_stayNRew',latencyType{nl}))) mean(eval(sprintf('latency_m.%s_switchNRew',latencyType{nl})))];
    
    sem = [std(eval(sprintf('latency_f.%s_stayNRew',latencyType{nl})))./sqrt(nf) std(eval(sprintf('latency_f.%s_switchNRew',latencyType{nl})))./sqrt(nf);...
        std(eval(sprintf('latency_m.%s_stayNRew',latencyType{nl}))./sqrt(nm)) std(eval(sprintf('latency_m.%s_switchNRew',latencyType{nl})))./sqrt(nm)];
    
    subplot(1,3,3); hold on
    b=bar(mu');
     b(1).EdgeColor = plotParams.femaleC;
    b(1).FaceColor = 'none';
    b(2).EdgeColor = plotParams.maleC;
    b(2).FaceColor = 'none';
    b(2).LineWidth = 1.5;
    b(1).LineWidth = 1.5;
    pause(.001);
    errorbar(b(1).XData+b(1).XOffset,mu(1,:),sem(1,:), 'Color',plotParams.femaleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    errorbar(b(2).XData+b(2).XOffset,mu(2,:),sem(2,:), 'Color',plotParams.maleC,'LineStyle','none','CapSize',0,'LineWidth',1.5)
    scatter(repmat(b(1).XData(1)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_stayNRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(1).XData(2)+b(1).XOffset,1,nf)-datasample(.02:.01:.1,nf), eval(sprintf('latency_f.%s_switchNRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(1)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_stayNRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    scatter(repmat(b(2).XData(2)+b(2).XOffset,1,nm)-datasample(.02:.01:.1,nm), eval(sprintf('latency_m.%s_switchNRew',latencyType{nl})), 10, 'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)

    set(gca,'XTick',b(1).XData,'XTickLabel',{'stay';'switch'},'YLim',[0 15])
    ylabel(latencyType{nl})
    title('Previously unrewarded trials')
    
    
    
    %% Stats
    
   
     thisLatency = cat(1,cell2mat(eval(sprintf('latency_f.%s_stayNRew_trials',latencyType{nl}))'),...
        cell2mat(eval(sprintf('latency_f.%s_switchNRew_trials',latencyType{nl}))'),...
        cell2mat(eval(sprintf('latency_m.%s_stayNRew_trials',latencyType{nl}))'),...
        cell2mat(eval(sprintf('latency_m.%s_switchNRew_trials',latencyType{nl}))'));
    Sex = cat(1,ones(size(cell2mat(eval(sprintf('latency_f.%s_stayNRew_trials',latencyType{nl}))'))),...
        ones(size(cell2mat(eval(sprintf('latency_f.%s_switchNRew_trials',latencyType{nl}))'))),...
        zeros(size(cell2mat(eval(sprintf('latency_m.%s_stayNRew_trials',latencyType{nl}))'))),...
        zeros(size(cell2mat(eval(sprintf('latency_m.%s_switchNRew_trials',latencyType{nl}))'))));
    Stay = cat(1,ones(size(cell2mat(eval(sprintf('latency_f.%s_stayNRew_trials',latencyType{nl}))'))),...
        zeros(size(cell2mat(eval(sprintf('latency_f.%s_switchNRew_trials',latencyType{nl}))'))),...
        ones(size(cell2mat(eval(sprintf('latency_m.%s_stayNRew_trials',latencyType{nl}))'))),...
        zeros(size(cell2mat(eval(sprintf('latency_m.%s_switchNRew_trials',latencyType{nl}))'))));
    
    AnimalID = cellfun(@(x,y) y.*ones(size(x))', eval(sprintf('latency_f.%s_stayNRew_trials',latencyType{nl})),num2cell(1:nf),'UniformOutput',false);
    AnimalID = cat(2,AnimalID, cellfun(@(x,y) y.*ones(size(x))', eval(sprintf('latency_f.%s_switchNRew_trials',latencyType{nl})),num2cell(1:nf),'UniformOutput',false));
    AnimalID = cat(2,AnimalID, cellfun(@(x,y) y.*ones(size(x))', eval(sprintf('latency_m.%s_stayNRew_trials',latencyType{nl})),num2cell(1:nm),'UniformOutput',false));
    AnimalID = cat(2,AnimalID, cellfun(@(x,y) y.*ones(size(x))', eval(sprintf('latency_m.%s_switchNRew_trials',latencyType{nl})),num2cell(1:nm),'UniformOutput',false));
    
    AnimalID = cell2mat(AnimalID)';
    
    X = table((thisLatency),'VariableNames',{'Latency'});
    X.Sex = nominal(Sex);
    X.Stay = nominal(Stay);
    X.AnimalID = nominal(AnimalID); 
    
    f = 'Latency ~ Sex*Stay + (1+Stay|AnimalID)';
    
    glme = fitlme(X,f,'DummyVarCoding','effects');
    
    stats.glme = fitlme(X,f,'DummyVarCoding','effects');
    stats.anova = dataset2table(anova(stats.glme,'DFMethod','Satterthwaite'));
    
    eval(sprintf('[stats.postHoc_%s_norew_switch_p,stats.postHoc_%s_norew_switch_h,stats.postHoc_%s_norew_switch_stats] = ranksum(latency_f.%s_switchNRew,latency_m.%s_switchNRew);',latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl})); 
    eval(sprintf('[stats.postHoc_%s_norew_stay_p,stats.postHoc_%s_norew_stay_h,stats.postHoc_%s_norew_stay_stats] = ranksum(latency_f.%s_stayNRew,latency_m.%s_stayNRew);',latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl})); 
    eval(sprintf('[stats.postHoc_%s_norew_female_p,stats.postHoc_%s_norew_female_h,stats.postHoc_%s_norew_female_stats] = signrank(latency_f.%s_switchNRew,latency_f.%s_stayNRew);',latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl})); 
    eval(sprintf('[stats.postHoc_%s_norew_male_p,stats.postHoc_%s_norew_male_h,stats.postHoc_%s_norew_male_stats] = signrank(latency_m.%s_switchNRew,latency_m.%s_stayNRew);',latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl},latencyType{nl})); 

end