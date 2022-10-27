function stats=latencyQuant_opto_old(cohort,ext,zscoreFlag,saveLoc,qFile,valType_plot,binNum,laserType)

% % Load data 
try
     load(fullfile(saveLoc, sprintf('%s_allData_opto_%s_zscore%d_bin%d_%s',cohort,ext,zscoreFlag,binNum,qFile)));
catch
    extractData_opto(cohort,ext,qFile,zscoreFlag,binNum,saveLoc);
    load(fullfile(saveLoc, sprintf('%s_allData_opto_%s_zscore%d_bin%d_%s',cohort,ext,zscoreFlag,binNum,qFile)));
end

ids_m = unique(Animal(Opsin==1&Sex==0));
ids_f = unique(Animal(Opsin==1&Sex==1));
ids_f_yfp = unique(Animal(Opsin==0&Sex==1));
ids_m_yfp = unique(Animal(Opsin == 0 &Sex==0)); 
ids = cat(1,ids_m,ids_f,ids_f_yfp,ids_m_yfp);


Latency(Trial==1) = NaN; % remove first trial of the session 


Latency_thresh = Latency;

for na = 1:numel(ids)
    thisSession = unique(Session(strcmp(Animal,ids{na})));
    for ns = 1:numel(thisSession)
        sessionIdx = find(Session==thisSession(ns)&strcmp(Animal,ids{na})); 
        threshIdx = sessionIdx(find(Latency(sessionIdx)>intervalThresh,1,'first'));
        Latency_thresh(threshIdx:sessionIdx(end)) = NaN;
    end
end
Latency = Latency_thresh;
% 
% idx = find(Laser==1);
% Laser(idx-1) = -1;
% Laser(idx) = 1;
% idx = find(Laser==2);
% Laser(idx-1) = -2;
% idx = find(Laser==3);
% Laser(idx+1) = -3; 

% 
% Laser_outcome = nan(size(Laser));
% Laser_outcome(Laser==1) = 1;
% Laser_outcome(Laser==0) = 0;
% Laser_np = nan(size(Laser));
% Laser_np(Laser==2) = 1;
% Laser_np(Laser==0) = 0; 
% Laser_ITI = nan(size(Laser));
% Laser_ITI(Laser == 3)=1; 
% Laser_ITI(Laser==0) = 0;

thisValue_ptile = Value_ptile{ strcmp(valType,valType_plot)}; % which type of value
thisValue = Value{ strcmp(valType,valType_plot)}; % which type of value
thisValue_zscore = Value_zscore{ strcmp(valType,valType_plot)}; % which type of value


Latency(Latency<.01) = 0.01; % correct for different temporal resolutions

%% Extract animal means

groups = {'f';'m';'f_yfp';'m_yfp'};

for ng = 1:numel(groups)
    thisIDs = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIDs)
        for nb = 1:binNum
            eval(sprintf('X_%s(na,nb) = nanmean(Latency(thisValue_ptile==nb&Animal==thisIDs(na)&Laser==0&SessionType==1));',groups{ng}))
            eval(sprintf('X_laser_%s(na,nb) = nanmean(Latency(thisValue_ptile==nb&Animal==thisIDs(na)&Laser==laserType&SessionType==1));',groups{ng}))
        end
    end
end


%% Female laser v no laser by value
for nb = 1:binNum
    %if lillietest(X_f(:,nb)) | lillietest(X_laser_f(:,nb))
        whichTest{1}(nb) = 1;
        [stats.pvals_posthoc(nb,1),~,stats.posthoc{nb,1}] = signrank(X_f(:,nb), X_laser_f(:,nb));
%     else
%         [~,stats.pvals_posthoc(nb,1),stats.posthoc{1,nb}] = ttest(X_f(:,nb), X_laser_f(:,nb));
%         whichTest{1}(nb) = 0; 
%     end
end

% Male laser v no laser by value
for nb = 1:binNum
   % if lillietest(X_m(:,nb)) | lillietest(X_laser_m(:,nb))
        whichTest{2}(nb) = 1;
        [stats.pvals_posthoc(nb,2),~,stats.posthoc{nb,2}] = signrank(X_m(:,nb), X_laser_m(:,nb));
%     else
%         [~,pvals(nb,2),stats{2,nb}] = ttest(X_m(:,nb), X_laser_m(:,nb));
%         whichTest{2}(nb) = 0; 
%     end
end

% female yfp
for nb = 1:binNum
    %if lillietest(X_f(:,nb)) | lillietest(X_laser_f(:,nb))
        whichTest{3}(nb) = 1;
        [stats.pvals_posthoc(nb,3),~,stats.posthoc{nb,3}] = signrank(X_f_yfp(:,nb), X_laser_f_yfp(:,nb) );
%     else
%         [~,stats.pvals_posthoc(nb,1),stats.posthoc{1,nb}] = ttest(X_f(:,nb), X_laser_f(:,nb));
%         whichTest{1}(nb) = 0; 
%     end
end

% Male laser v no laser by value
for nb = 1:binNum
   % if lillietest(X_m(:,nb)) | lillietest(X_laser_m(:,nb))
        whichTest{4}(nb) = 1;
        [stats.pvals_posthoc(nb,4),~,stats.posthoc{nb,4}] = signrank(X_m_yfp(:,nb), X_laser_m_yfp(:,nb));
%     else
%         [~,pvals(nb,2),stats{2,nb}] = ttest(X_m(:,nb), X_laser_m(:,nb));
%         whichTest{2}(nb) = 0; 
%     end
end

%% Mixed-effects regession of  difference between laser and no laser

value = [];
LatencyDiff = [];
sex = [];
opsin = [];
animal = [];
for nb = 1:numel(unique(thisValue_ptile))
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




%% Plot 
    
    
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;

   
  


figure();
% Plot female nphr outcome laser v no laser
subplot(2,2,1); hold on
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_f,1));%+datasample(0:.01:.25,size(X_f,1));
   scatter(xaxis, X_f(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_f,1));%-datasample(0:.01:.25,size(X_f,1));
   scatter(xaxis, X_laser_f(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
p=errorbar(nanmean(X_laser_f),nanstd(X_laser_f)./sqrt(size(X_laser_f,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(1:numel(unique(thisValue_ptile)),nanmean(X_f),nanstd(X_f)./sqrt(size(X_f,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
legend(p,{'laser';'no laser'})

subplot(2,2,2); hold on
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_m,1));%+datasample(0:.01:.25,size(X_m,1));
   scatter(xaxis, X_m(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_m,1));%-datasample(0:.01:.25,size(X_m,1));
   scatter(xaxis, X_laser_m(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
p = errorbar(nanmean(X_m),nanstd(X_m)./sqrt(size(X_m,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0);
p(2)=errorbar(nanmean(X_laser_m),nanstd(X_laser_m)./sqrt(size(X_laser_m,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0);
legend(p,{'laser';'no laser'})

subplot(2,2,3); hold on
errorbar(nanmean(X_f_yfp),nanstd(X_f_yfp)./sqrt(size(X_f_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(X_laser_f_yfp),nanstd(X_laser_f_yfp)./sqrt(size(X_laser_f_yfp,1)), 'Color',femaleC,'LineWidth',1.5,'CapSize',0)
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_f_yfp,1));%+datasample(0:.01:.25,size(X_f_yfp,1));
   scatter(xaxis, X_f_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_f_yfp,1));%-datasample(0:.01:.25,size(X_f_yfp,1));
   scatter(xaxis, X_laser_f_yfp(:,nb), 15, 'MarkerFaceColor', femaleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end

subplot(2,2,4); hold on
errorbar(nanmean(X_m_yfp),nanstd(X_m_yfp)./sqrt(size(X_m_yfp,1)), 'Color',[.3 .3 .3],'LineWidth',1.5,'CapSize',0)
errorbar(nanmean(X_laser_m_yfp),nanstd(X_laser_m_yfp)./sqrt(size(X_laser_m_yfp,1)), 'Color',maleC,'LineWidth',1.5,'CapSize',0)
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_m_yfp,1));%+datasample(0:.01:.25,size(X_m_yfp,1));
   scatter(xaxis, X_m_yfp(:,nb), 15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end
for nb = 1:numel(unique(thisValue_ptile))
   xaxis = repmat(nb,1,size(X_m_yfp,1));%-datasample(0:.01:.25,size(X_m_yfp,1));
   scatter(xaxis, X_laser_m_yfp(:,nb), 15, 'MarkerFaceColor', maleC, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',.5);
end


for np = 1:4
    subplot(2,2,np)
    axis square
    ylabel('Trial initiation latency (sec)')
    set(gca,'YLim',[0,20],'XLim', [.5 max(unique(thisValue_ptile))+.5])
    if laserType==3
        set(gca,'YLim',[0 20],'XLim', [.5 max(unique(thisValue_ptile))+.5])
    end
    xlabel(sprintf('%s percentile',valType_plot))
end



