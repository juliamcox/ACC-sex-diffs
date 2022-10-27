function plotPsychCurves(sessionLength, perfThresh, qFile, binNum, aids, femaleFlag)

% Plot real and estimated choice x value plots for each animal 

%% Parameters
qvals = linspace(-1,1,binNum+1); % q value bins for estimating choice
stdInterval = normcdf(1, 0, 1) - normcdf(-1, 0, 1);
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); 
if femaleFlag
    plotC = plotParams.femaleC;
else
    plotC = plotParams.maleC; 
end

numSubs = [4, 4];
panelCounter = 1;
figCounter(1) = 1;
figCounter(2) = 2; 

savename = fullfile(whereAreWe('figurecode'),'processed_data','psychCurves',sprintf('qDiff_choice_%s.ps',date));
savename2 = fullfile(whereAreWe('figurecode'),'processed_data','psychCurves',sprintf('qDiffQuant_choice_%s.ps',date));
%% Estimate choice by session and run for na = 1:numel(aids)
for na = 1:numel(aids)
    try
        load(fullfile(whereAreWe('bucket'),'Operant',aids{na},sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile))); % load q-values and choice
    catch
        valueExtraction_TAB(aids(na),qFile,sessionLength,perfThresh);
        load(fullfile(whereAreWe('bucket'),'Operant',aids{na},sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile))); % load q-values and choice
   end
    
    load(fullfile(whereAreWe('bucket'),'Operant',aids{na},qFile));
    
    sessionIDs = find(contains({qLearn.fList(:).name},flist));
    idx = boolean(sum(qLearn.session==sessionIDs,2));
    qDiff = qLearn.QDiff(idx); 
    predChoice = qLearn.estChoice(:,idx); 
    realChoice = qLearn.choice(idx)==0; 
    
    
    
    % bin trials by q right - q left (qDiff)
    [~,~,qDiff_bin] = histcounts(qDiff,qvals);
    [~,~,qDiffQuant_bin] = histcounts(qDiff,prctile(qDiff,linspace(0,100,10)));
    
    try
        for nb = 1:numel(qvals)-1
            qDiff_predChoice(na,nb) = sum(predChoice(qDiff_bin==nb))./sum(qDiff_bin==nb);
            [~,temp] = binofit(sum(predChoice(qDiff_bin==nb)), sum(qDiff_bin==nb));
            %             [~, temp] = binointerval(sum(predChoice(qDiff_bin==nb)), sum(qDiff_bin==nb), 1-stdInterval);
            qDiff_predChoice_lower(na,nb) = temp(1);
            qDiff_predChoice_upper(na,nb) = temp(2);
            
            qDiff_realChoice(na,nb) = sum(realChoice(qDiff_bin==nb))./sum(qDiff_bin==nb);
            [~,temp] = binofit(sum(realChoice(qDiff_bin==nb)), sum(qDiff_bin==nb));
            %             [~, temp] = binointerval(sum(realChoice(qDiff_bin==nb)), sum(qDiff_bin==nb), 1-stdInterval);
            qDiff_realChoice_lower(na,nb) = temp(1);
            qDiff_realChoice_upper(na,nb) = temp(2);
            %
            qDiffQuant_predChoice(na,nb) = sum(predChoice(qDiffQuant_bin==nb))./sum(qDiffQuant_bin==nb);
            [~,temp] = binofit(sum(predChoice(qDiffQuant_bin==nb)), sum(qDiffQuant_bin==nb));
            %             [~, temp] = binointerval(sum(predChoice(qDiffQuant_bin==nb)), sum(qDiffQuant_bin==nb), 1-stdInterval);
            qDiffQuant_predChoice_lower(na,nb) = temp(1);
            qDiffQuant_predChoice_upper(na,nb) = temp(2);
            
            
            qDiffQuant_realChoice(na,nb) = sum(realChoice(qDiffQuant_bin==nb))./sum(qDiffQuant_bin==nb);
            [~,temp] = binofit(sum(realChoice(qDiffQuant_bin==nb)), sum(qDiffQuant_bin==nb));
            %             [~, temp] = binointerval(sum(realChoice(qDiffQuant_bin==nb)), sum(qDiffQuant_bin==nb), 1-stdInterval);
            qDiffQuant_realChoice_lower(na,nb) = temp(1);
            qDiffQuant_realChoice_upper(na,nb) = temp(2);
            
        end
    catch
        keyboard
    end
    
    if panelCounter <= numSubs(1)*numSubs(2)
        figure(figCounter(1))
        subplot(numSubs(1),numSubs(2),panelCounter); hold on
        p(1)=plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice(na,:),'Color', plotC,'LineWidth',1.5);
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice_upper(na,:),'Color', [plotC .5],'LineWidth',1);
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice_lower(na,:),'Color', [plotC .5],'LineWidth',1);

        p(2)=plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_predChoice(na,:),'Color', [.7 .7 .7],'LineWidth',1.5);    
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_predChoice_upper(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        plot(qvals(1:end-1)+mean(diff(qvals))/2,qDiff_predChoice_lower(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        title(aids{na})
        ylabel('Prob(R)')
        xlabel('QR-QL')
        if panelCounter==1
            [l,icons]=legend(p,{'real';'estimated'},'Location','SouthEast','Box','off');
            icons(3).XData = [0.3 0.4198];
            icons(5).XData = [0.3 0.4198];
            l.Position(1) = l.Position(1)*1.25;
        end
        set(gca,'YLim',[0 1])

        figure(figCounter(2))
        subplot(numSubs(1),numSubs(2),panelCounter); hold on
        p(1)=plot(1:numel(qvals)-1, qDiffQuant_realChoice(na,:),'Color', plotC,'LineWidth',1.5);   
        plot(1:numel(qvals)-1, qDiffQuant_realChoice_upper(na,:),'Color', [plotC .5],'LineWidth',1);
        plot(1:numel(qvals)-1, qDiffQuant_realChoice_lower(na,:),'Color', [plotC .5],'LineWidth',1);
        p(2)=plot(1:numel(qvals)-1, qDiffQuant_predChoice(na,:),'Color', [.7 .7 .7],'LineWidth',1.5);     
        plot(1:numel(qvals)-1, qDiffQuant_predChoice_upper(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        plot(1:numel(qvals)-1,qDiffQuant_predChoice_lower(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        title(aids{na})
        if panelCounter==1
            [l,icons]=legend(p,{'real';'estimated'},'Location','SouthEast','Box','off');
            icons(3).XData = [0.3 0.4198];
            icons(5).XData = [0.3 0.4198];
            l.Position(1) = l.Position(1)*1.25;
        end
        set(gca,'YLim',[0 1])
        ylabel('Prob(R)')
        xlabel('QR-QL quantile')
        axis square
        panelCounter = panelCounter+1; 
    else
%         figure(1)
%         print(savename,'-dpsc','-append')
%         figure(2)
%         print(savename2,'-dpsc','-append')
%         panelCounter = 1;
%         close all
        panelCounter = 1;
        figCounter = figCounter+2;
        figure(figCounter(1))
        subplot(numSubs(1),numSubs(2),panelCounter); hold on
        p(1)=plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice(na,:),'Color', plotC,'LineWidth',1.5);
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice_upper(na,:),'Color', [plotC .5],'LineWidth',1);
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_realChoice_lower(na,:),'Color', [plotC .5],'LineWidth',1);

        p(2)=plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_predChoice(na,:),'Color', [.7 .7 .7],'LineWidth',1.5);    
        plot(qvals(1:end-1)+mean(diff(qvals))/2, qDiff_predChoice_upper(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        plot(qvals(1:end-1)+mean(diff(qvals))/2,qDiff_predChoice_lower(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        title(aids{na})
        ylabel('Prob(R)')
        xlabel('QR-QL')
        if panelCounter==1
            [l,icons]=legend(p,{'real';'estimated'},'Location','SouthEast','Box','off');
            icons(3).XData = [0.3 0.4198];
            icons(5).XData = [0.3 0.4198];
            l.Position(1) = l.Position(1)*1.25;
        end
        set(gca,'YLim',[0 1])

        figure(figCounter(2))
        subplot(numSubs(1),numSubs(2),panelCounter); hold on
        p(1)=plot(1:numel(qvals)-1, qDiffQuant_realChoice(na,:),'Color', plotC,'LineWidth',1.5);   
        plot(1:numel(qvals)-1, qDiffQuant_realChoice_upper(na,:),'Color', [plotC .5],'LineWidth',1);
        plot(1:numel(qvals)-1, qDiffQuant_realChoice_lower(na,:),'Color', [plotC .5],'LineWidth',1);
        p(2)=plot(1:numel(qvals)-1, qDiffQuant_predChoice(na,:),'Color', [.7 .7 .7],'LineWidth',1.5);
        plot(1:numel(qvals)-1, qDiffQuant_predChoice_upper(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        plot(1:numel(qvals)-1,qDiffQuant_predChoice_lower(na,:),'Color', [.7 .7 .7 .5],'LineWidth',1);
        title(aids{na})
        if panelCounter==1
            [l,icons]=legend(p,{'real';'estimated'},'Location','SouthEast','Box','off');
            icons(3).XData = [0.3 0.4198];
            icons(5).XData = [0.3 0.4198];
            l.Position(1) = l.Position(1)*1.25;
        end
        set(gca,'YLim',[0 1])
        ylabel('Prob(R)')
        xlabel('QR-QL quantile')
        axis square
        panelCounter = panelCounter+1;
    end
end

figure(1);
for nf = 1:2:figCounter(1)
    figure(nf)
    print(savename,'-dpsc','-append')
end
for nf = 2:2:figCounter(2)
    figure(nf)
    print(savename2,'-dpsc','-append')
end
