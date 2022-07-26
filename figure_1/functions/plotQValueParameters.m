function stats=plotQValueParameters(params_f,params_m,qFile)

%%% Figure 1

% Calculate and plot q-value model parameter estimates for males and females (control sessions)

% Dependencies:
% requires mat files converted from mpc
% requires file with trial-by-trial Q-value estimates 
% functions: generateAnimalLists (to extract animal IDs to plot); whereAreWe (for generating filenames); getIndices (function to extract indices for laser, prev reward etc)


%% Parameters


% saveLoc       = fullfile(whereAreWe('bucket'),'Manuscript_figures', 'fig1', 'QLearn_params', sprintf('perfThresh%s_%s',num2str(perfThresh),ext));
% if ~isdir(saveLoc)
%     mkdir(saveLoc)
% end

a             = .01; % alpha for significance tests 

plotParams    = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); % load plot parameters 

%% Extract parameters for males and females



paramsTable = table(cat(1,params_f.alpha',params_m.alpha'),'VariableNames',{'alpha'});
paramsTable.alpha_sem = cat(1,params_f.alpha_sem',params_m.alpha_sem');
paramsTable.beta = (cat(1,params_f.beta',params_m.beta'));
paramsTable.beta_sem = cat(1,params_f.beta_sem',params_m.beta_sem');
paramsTable.stay = (cat(1,params_f.stay',params_m.stay'));
paramsTable.stay_sem = cat(1,params_f.stay_sem',params_m.stay_sem');
paramsTable.side = (cat(1,params_f.side',params_m.side'));
paramsTable.side_sem = cat(1,params_f.side_sem',params_m.side_sem');


paramsTable2 = table(arrayfun(@(x,y) sprintf('%.2f \x00B1 %.2f',x,y),cat(1,params_f.alpha',params_m.alpha'),cat(1,params_f.alpha_sem',params_m.alpha_sem'),'UniformOutput',false),'VariableNames',{'alpha'});
paramsTable2.beta = arrayfun(@(x,y) sprintf('%.2f \x00B1 %.2f',x,y),cat(1,params_f.beta',params_m.beta'),cat(1,params_f.beta_sem',params_m.beta_sem'),'UniformOutput',false);
paramsTable2.stay = arrayfun(@(x,y) sprintf('%.2f \x00B1 %.2f',x,y),cat(1,params_f.stay',params_m.stay'),cat(1,params_f.stay_sem',params_m.stay_sem'),'UniformOutput',false);
paramsTable2.side = arrayfun(@(x,y) sprintf('%.2f \x00B1 %.2f',x,y),cat(1,params_f.side',params_m.side'),cat(1,params_f.side_sem',params_m.side_sem'),'UniformOutput',false);

%% Plot parameters for males and females
if contains(qFile,'2alpha')
    mu_f  = [mean(params_f.alpha) mean(params_f.alpha2) mean(params_f.beta) mean(params_f.stay) mean(params_f.side)];
    mu_m  = [mean(params_m.alpha) mean(params_m.alpha2) mean(params_m.beta) mean(params_m.stay) mean(params_m.side)];
    sem_f = [std(params_f.alpha)./sqrt(numel(params_f.alpha)) std(params_f.alpha2)./sqrt(numel(params_f.alpha2)) std(params_f.beta)./sqrt(numel(params_f.beta)) std(params_f.stay)./sqrt(numel(params_f.stay)) std(params_f.side)./sqrt(numel(params_f.side))];
    sem_m = [std(params_m.alpha)./sqrt(numel(params_m.alpha)) std(params_m.alpha2)./sqrt(numel(params_m.alpha2)) std(params_m.beta)./sqrt(numel(params_m.beta)) std(params_m.stay)./sqrt(numel(params_m.stay)) std(params_m.side)./sqrt(numel(params_m.side))];
    
    vec = cat(1,params_f.alpha',params_f.alpha2',params_m.alpha',params_m.alpha2');
    groupvec = cat(1, [ones(size(params_f.alpha')) ones(size(params_f.alpha'))],[ones(size(params_f.alpha2')).*2 ones(size(params_f.alpha2'))],[ones(size(params_m.alpha')) zeros(size(params_m.alpha'))],[ones(size(params_m.alpha2')).*2 zeros(size(params_m.alpha2'))]);
    [stats1.p_anova,stats1.t_anova,stats1.s_anova] = anovan(vec, groupvec,'varnames', {'outcome';'sex'}, 'display', 'off', 'model','full');
    figure(); multcompare(stats1.s_anova,'dimension',[1 2]);
    
    
    figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
    b = bar([mu_f' mu_m']);
    b(1).EdgeColor = plotParams.femaleC;
    b(1).FaceColor = 'none';
    b(1).LineWidth = 1.5;
    b(2).EdgeColor = plotParams.maleC;
    b(2).FaceColor = 'none';
    b(2).LineWidth = 1.5;
    pause(.001);
    errorbar(b(1).XData+b(1).XOffset,b(1).YData,sem_f,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.femaleC)
    errorbar(b(2).XData+b(2).XOffset,b(2).YData,sem_m,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.maleC)
    
    % Plot individual data points
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha))+ones(size(params_f.alpha)).*b(1).XData(1)+b(1).XOffset;
    scatter(xaxis, params_f.alpha, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.alpha,params_m.alpha,'alpha',a)
        text('*',b(1).XData(1),3)
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha2))+ones(size(params_f.alpha2)).*b(1).XData(2)+b(1).XOffset;
    scatter(xaxis, params_f.alpha2, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.alpha2,params_m.alpha2,'alpha',a)
        text('*',b(1).XData(1),3)
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.beta))+ones(size(params_f.beta)).*b(1).XData(3)+b(1).XOffset;
    scatter(xaxis, params_f.beta, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.beta,params_m.beta,'alpha',a)
        text('*',b(1).XData(2),3)
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.stay))+ones(size(params_f.stay)).*b(1).XData(4)+b(1).XOffset;
    scatter(xaxis, params_f.stay, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.stay,params_m.stay,'alpha',a)
        text(b(1).XData(3),3,'*')
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.side))+ones(size(params_f.side)).*b(1).XData(5)+b(1).XOffset;
    scatter(xaxis, params_f.side, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.side,params_m.side,'alpha',a)
        text('*',b(1).XData(4),3)
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha))+ones(size(params_m.alpha)).*b(2).XData(1)+b(2).XOffset;
    scatter(xaxis, params_m.alpha, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha2))+ones(size(params_m.alpha2)).*b(2).XData(2)+b(2).XOffset;
    scatter(xaxis, params_m.alpha2, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.beta))+ones(size(params_m.beta)).*b(2).XData(3)+b(2).XOffset;
    scatter(xaxis, params_m.beta, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.stay))+ones(size(params_m.stay)).*b(2).XData(4)+b(2).XOffset;
    scatter(xaxis, params_m.stay, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.side))+ones(size(params_m.side)).*b(2).XData(5)+b(2).XOffset;
    scatter(xaxis, params_m.side, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    
    
    set(gca,'XTick',b(1).XData, 'XTickLabel', {'alpha_pos';'alpha_neg';'beta';'stay';'bias'})
    ylabel('Parameter estimate')
    
else
    mu_f  = [mean(params_f.alpha) mean(params_f.beta) mean(params_f.stay) mean(params_f.side)];
    mu_m  = [mean(params_m.alpha) mean(params_m.beta) mean(params_m.stay) mean(params_m.side)];
    sem_f = [std(params_f.alpha)./sqrt(numel(params_f.alpha)) std(params_f.beta)./sqrt(numel(params_f.beta)) std(params_f.stay)./sqrt(numel(params_f.stay)) std(params_f.side)./sqrt(numel(params_f.side))];
    sem_m = [std(params_m.alpha)./sqrt(numel(params_m.alpha)) std(params_m.beta)./sqrt(numel(params_m.beta)) std(params_m.stay)./sqrt(numel(params_m.stay)) std(params_m.side)./sqrt(numel(params_m.side))];
    
    figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
    b = bar([mu_f' mu_m']);
    b(1).EdgeColor = plotParams.femaleC;
    b(1).FaceColor = 'none';
    b(1).LineWidth = 1.5;
    b(2).EdgeColor = plotParams.maleC;
    b(2).FaceColor = 'none';
    b(2).LineWidth = 1.5;
    pause(.001);
    % Plot errorbars
    errorbar(b(1).XData+b(1).XOffset,b(1).YData,sem_f,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.femaleC)
    errorbar(b(2).XData+b(2).XOffset,b(2).YData,sem_m,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.maleC)
    
    % Plot individual data points
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha))+ones(size(params_f.alpha)).*b(1).XData(1)+b(1).XOffset;
    scatter(xaxis, params_f.alpha, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.alpha,params_m.alpha,'alpha',a)
        text(b(1).XData(1),3,'*')
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.beta))+ones(size(params_f.beta)).*b(1).XData(2)+b(1).XOffset;
    scatter(xaxis, params_f.beta, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.beta,params_m.beta,'alpha',a)
        text(b(1).XData(2),3,'*')
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.stay))+ones(size(params_f.stay)).*b(1).XData(3)+b(1).XOffset;
    scatter(xaxis, params_f.stay, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.stay,params_m.stay,'alpha',a)
        text(b(1).XData(3),3,'*')
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.side))+ones(size(params_f.side)).*b(1).XData(4)+b(1).XOffset;
    scatter(xaxis, params_f.side, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    if ttest2(params_f.side,params_m.side,'alpha',a)
        text(b(1).XData(4),3,'*')
    end
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha))+ones(size(params_m.alpha)).*b(2).XData(1)+b(2).XOffset;
    scatter(xaxis, params_m.alpha, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.beta))+ones(size(params_m.beta)).*b(2).XData(2)+b(2).XOffset;
    scatter(xaxis, params_m.beta, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.stay))+ones(size(params_m.stay)).*b(2).XData(3)+b(2).XOffset;
    scatter(xaxis, params_m.stay, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.side))+ones(size(params_m.side)).*b(2).XData(4)+b(2).XOffset;
    scatter(xaxis, params_m.side, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
    
    [~,stats.ttest_alpha_p,~,stats.ttest_alpha_stat] = ttest2(params_f.alpha,params_m.alpha);
    [~,stats.ttest_stay_p,~,stats.ttest_stay_stat] = ttest2(params_f.stay,params_m.stay);
    [~,stats.ttest_beta_p,~,stats.ttest_beta_stat] = ttest2(params_f.beta,params_m.beta);
    [~,stats.ttest_side_p,~,stats.ttest_side_stat] = ttest2(params_f.side,params_m.side);

    
    set(gca,'XTick', b(1).XData, 'XTickLabel', {'alpha';'beta';'stay';'bias'})
    ylabel('Parameter estimate','FontSize',10)
    
    figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
    female = cat(1,ones(size(params_f.alpha')),2.*ones(size(params_m.alpha')),3.*ones(size(params_f.alpha')),4.*ones(size(params_m.alpha')),5.*ones(size(params_f.alpha')),6.*ones(size(params_m.alpha')),7.*ones(size(params_f.alpha')),8.*ones(size(params_m.alpha')));
    b=boxplot(cat(1,params_f.alpha',params_m.alpha',params_f.beta',params_m.beta',params_f.stay',params_m.stay',params_f.side',params_m.side'),...
        female,'Labels',{'\alpha';'{alpha}';'beta';'beta';'stay';'stay';'side';'side'},'Color',cat(1,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
    ax = gca;
    ax.TickLabelInterpreter = 'tex';

    ax.XTickLabel= {'\alpha';[];'\beta_{value}';[];'\beta_{stay}';[];'\beta_{side}';[]};

    box off
    
    
end
%print(fullfile(saveLoc,'qParams_session.pdf'),'-dpdf');

%% Make table 
tbl_f = table([1:numel(params_f.alpha)]',repmat({'female'},numel(params_f.alpha),1),'VariableNames',{'Animal';'Sex'});
tbl_f.alpha = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.alpha,params_f.alpha_sem,'UniformOutput',false);
tbl_f.beta_value = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.beta,params_f.beta_sem,'UniformOutput',false)';
tbl_f.beta_stay = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.stay,params_f.stay_sem,'UniformOutput',false);
tbl_f.beta_side = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.side,params_f.side_sem,'UniformOutput',false);

tbl_m = table([numel(params_f.alpha)+1:numel(params_m.alpha)+numel(params_f.alpha)]',repmat({'female'},numel(params_m.alpha),1),'VariableNames',{'Animal';'Sex'});
tbl_m.alpha = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.3f')), params_m.alpha,params_m.alpha_sem,'UniformOutput',false);
tbl_m.beta_value = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.beta,params_m.beta_sem,'UniformOutput',false)';
tbl_m.beta_stay = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.stay,params_m.stay_sem,'UniformOutput',false)';
tbl_m.beta_side = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.side,params_m.side_sem,'UniformOutput',false)';

writetable(tbl_f,fullfile(whereAreWe('figurecode'),'processed_data','stats_tables','qmodel_f.csv'))
writetable(tbl_m,fullfile(whereAreWe('figurecode'),'processed_data','stats_tables','qmodel_m.csv'))


%% Plot parameters for males and females (animal level)
% if contains(qFile,'2alpha')
%     mu_f  = [mean(params_f.alpha_m) mean(params_f.alpha2_m) mean(params_f.beta_m) mean(params_f.stay_m) mean(params_f.side_m)];
%     mu_m  = [mean(params_m.alpha_m) mean(params_m.alpha2_m) mean(params_m.beta_m) mean(params_m.stay_m) mean(params_m.side_m)];
%     sem_f = [std(params_f.alpha_m)./sqrt(numel(params_f.alpha_m)) std(params_f.alpha2_m)./sqrt(numel(params_f.alpha2_m)) std(params_f.beta_m)./sqrt(numel(params_f.beta_m)) std(params_f.stay_m)./sqrt(numel(params_f.stay_m)) std(params_f.side_m)./sqrt(numel(params_f.side_m))];
%     sem_m = [std(params_m.alpha_m)./sqrt(numel(params_m.alpha_m)) std(params_m.alpha2_m)./sqrt(numel(params_m.alpha2_m)) std(params_m.beta_m)./sqrt(numel(params_m.beta_m)) std(params_m.stay_m)./sqrt(numel(params_m.stay_m)) std(params_m.side_m)./sqrt(numel(params_m.side_m))];
%     
%     vec = cat(1,params_f.alpha_m',params_f.alpha2_m',params_m.alpha_m',params_m.alpha2_m');
%     groupvec = cat(1, [ones(size(params_f.alpha_m')) ones(size(params_f.alpha_m'))],[ones(size(params_f.alpha2_m')).*2 ones(size(params_f.alpha2_m'))],[ones(size(params_m.alpha_m')) zeros(size(params_m.alpha_m'))],[ones(size(params_m.alpha2_m')).*2 zeros(size(params_m.alpha2_m'))]);
%     [stats.p_anova,stats.t_anova,stats.s_anova] = anovan(vec, groupvec,'varnames', {'outcome';'sex'}, 'display', 'off', 'model','full');
%     figure(); multcompare(stats.s_anova,'dimension',[1 2]);
%     
%     figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
%     b = bar([mu_f' mu_m']);
%     b(1).EdgeColor = plotParams.femaleC;
%     b(1).FaceColor = 'none';
%     b(1).LineWidth = 1.5;
%     b(2).EdgeColor = plotParams.maleC;
%     b(2).FaceColor = 'none';
%     b(2).LineWidth = 1.5;
%     pause(.001);
%     errorbar(b(1).XData+b(1).XOffset,b(1).YData,sem_f,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.femaleC)
%     errorbar(b(2).XData+b(2).XOffset,b(2).YData,sem_m,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.maleC)
%     
%     % Plot individual data points
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha_m))+ones(size(params_f.alpha_m)).*b(1).XData(1)+b(1).XOffset;
%     scatter(xaxis, params_f.alpha_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.alpha_m,params_m.alpha_m,'alpha',a)
%         text('*',b(1).XData(1),3)
%     end
%     
%      xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha2_m))+ones(size(params_f.alpha2_m)).*b(1).XData(2)+b(1).XOffset;
%     scatter(xaxis, params_f.alpha2_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.alpha2_m,params_m.alpha2_m,'alpha',a)
%         text('*',b(1).XData(1),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.beta_m))+ones(size(params_f.beta_m)).*b(1).XData(3)+b(1).XOffset;
%     scatter(xaxis, params_f.beta_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.beta_m,params_m.beta_m,'alpha',a)
%         text('*',b(1).XData(2),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.stay_m))+ones(size(params_f.stay_m)).*b(1).XData(4)+b(1).XOffset;
%     scatter(xaxis, params_f.stay_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.stay_m,params_m.stay_m,'alpha',a)
%         text(b(1).XData(3),3,'*')
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.side_m))+ones(size(params_f.side_m)).*b(1).XData(5)+b(1).XOffset;
%     scatter(xaxis, params_f.side_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.side_m,params_m.side_m,'alpha',a)
%         text('*',b(1).XData(4),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha_m))+ones(size(params_m.alpha_m)).*b(2).XData(1)+b(2).XOffset;
%     scatter(xaxis, params_m.alpha_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha2_m))+ones(size(params_m.alpha2_m)).*b(2).XData(2)+b(2).XOffset;
%     scatter(xaxis, params_m.alpha2_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.beta_m))+ones(size(params_m.beta_m)).*b(2).XData(3)+b(2).XOffset;
%     scatter(xaxis, params_m.beta_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.stay_m))+ones(size(params_m.stay_m)).*b(2).XData(4)+b(2).XOffset;
%     scatter(xaxis, params_m.stay_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.side_m))+ones(size(params_m.side_m)).*b(2).XData(5)+b(2).XOffset;
%     scatter(xaxis, params_m.side_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     
%     
%   
%     
%     set(gca,'XTick',b(1).XData, 'XTickLabel', {'alpha_pos';'alpha_neg';'beta';'stay';'bias'})
%     ylabel('Parameter estimate')
%     
%     print(fullfile(saveLoc,'qParams_animal.pdf'),'-dpdf');
% else
%     mu_f  = [mean(params_f.alpha_m) mean(params_f.beta_m) mean(params_f.stay_m) mean(params_f.side_m)];
%     mu_m  = [mean(params_m.alpha_m) mean(params_m.beta_m) mean(params_m.stay_m) mean(params_m.side_m)];
%     sem_f = [std(params_f.alpha_m)./sqrt(numel(params_f.alpha_m)) std(params_f.beta_m)./sqrt(numel(params_f.beta_m)) std(params_f.stay_m)./sqrt(numel(params_f.stay_m)) std(params_f.side_m)./sqrt(numel(params_f.side_m))];
%     sem_m = [std(params_m.alpha_m)./sqrt(numel(params_m.alpha_m)) std(params_m.beta_m)./sqrt(numel(params_m.beta_m)) std(params_m.stay_m)./sqrt(numel(params_m.stay_m)) std(params_m.side_m)./sqrt(numel(params_m.side_m))];
%     
%     figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
%     b = bar([mu_f' mu_m']);
%     b(1).EdgeColor = plotParams.femaleC;
%     b(1).FaceColor = 'none';
%     b(1).LineWidth = 1.5;
%     b(2).EdgeColor = plotParams.maleC;
%     b(2).FaceColor = 'none';
%     b(2).LineWidth = 1.5;
%     pause(.001);
%     errorbar(b(1).XData+b(1).XOffset,b(1).YData,sem_f,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.femaleC)
%     errorbar(b(2).XData+b(2).XOffset,b(2).YData,sem_m,'LineStyle', 'none','CapSize',0,'LineWidth',1.5,'Color',plotParams.maleC)
%     
%     % Plot individual data points
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.alpha_m))+ones(size(params_f.alpha_m)).*b(1).XData(1)+b(1).XOffset;
%     scatter(xaxis, params_f.alpha_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.alpha_m,params_m.alpha_m,'alpha',a)
%         text('*',b(1).XData(1),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.beta_m))+ones(size(params_f.beta_m)).*b(1).XData(2)+b(1).XOffset;
%     scatter(xaxis, params_f.beta_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.beta_m,params_m.beta_m,'alpha',a)
%         text('*',b(1).XData(2),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.stay_m))+ones(size(params_f.stay_m)).*b(1).XData(3)+b(1).XOffset;
%     scatter(xaxis, params_f.stay_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.stay_m,params_m.stay_m,'alpha',a)
%         text(b(1).XData(3),3,'*')
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_f.side_m))+ones(size(params_f.side_m)).*b(1).XData(4)+b(1).XOffset;
%     scatter(xaxis, params_f.side_m, 8,'MarkerFaceColor',[plotParams.femaleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     if ttest2(params_f.side_m,params_m.side_m,'alpha',a)
%         text('*',b(1).XData(4),3)
%     end
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.alpha_m))+ones(size(params_m.alpha_m)).*b(2).XData(1)+b(2).XOffset;
%     scatter(xaxis, params_m.alpha_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.beta_m))+ones(size(params_m.beta_m)).*b(2).XData(2)+b(2).XOffset;
%     scatter(xaxis, params_m.beta_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.stay_m))+ones(size(params_m.stay_m)).*b(2).XData(3)+b(2).XOffset;
%     scatter(xaxis, params_m.stay_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     xaxis = -.05+datasample(-.025:.001:.025,numel(params_m.side_m))+ones(size(params_m.side_m)).*b(2).XData(4)+b(2).XOffset;
%     scatter(xaxis, params_m.side_m, 8,'MarkerFaceColor',[plotParams.maleC],'MarkerEdgeColor','none','MarkerFaceAlpha',plotParams.barAlpha)
%     
%     if contains(qFile,'2alpha')
%     end
%     
%     set(gca,'XTick',b(1).XData, 'XTickLabel', {'alpha';'beta';'stay';'bias'})
%     ylabel('Parameter estimate')
%     
%    % print(fullfile(saveLoc,'qParams_animal.pdf'),'-dpdf');
% end
end


