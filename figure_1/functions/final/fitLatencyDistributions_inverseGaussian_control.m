function [stats,fits_f,fits_m] = fitLatencyDistributions_inverseGaussian_control(behaviorTable,latencyType)

%behaviorTable = trial by trial data
%latencyType = which latency to fit (trialInit, trialInit_thresh,leverPress

%% Parameters

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

behaviorTable = behaviorTable(behaviorTable.laserSession==0,:);
if contains(latencyType,'thresh')
   behaviorTable = behaviorTable(~isnan(behaviorTable.trialInit_thresh),:); 
end
Latency = eval(sprintf('behaviorTable.%s;',latencyType)); 
Latency(Latency<.01) = 0.01;


options.MaxIter = 2000; %increase number of iterations
options.MaxFunEvals = 2000; 

% get indices for each animal 
ids_f = unique(behaviorTable.aID(behaviorTable.female==1));
ids_m = unique(behaviorTable.aID(behaviorTable.female==0));


minVal = 2;
maxVal = 3;

numParams = 3;
thisParams = {'mu';'theta';'lambda'};



%% Fit shifted inverse gaussian pdf

custpdf = @(data,mu,shift,lambda) pdf('InverseGaussian',data-shift,mu,lambda);
custcdf = @(data,mu,shift,lambda) cdf('InverseGaussian',data-shift,mu,lambda);

groups = {'f';'m'};

for ng = 1:numel(groups)
    clear thisFits
    thisIds = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIds)
        
        data = Latency(strcmp(behaviorTable.aID,thisIds{na}));
        try
            thisFits{na} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[2 -.1 1],'lowerbound',[0 -2 0],'options',options);
        catch
            thisFits{na} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -2 0],'options',options);
        end
    end
    eval(sprintf('fits_%s = thisFits;',groups{ng}))
    clear thisFits
end


%% Plot individual fits
bins = [0:.3:10.3 300];
bins2 = [0:.3:10.5];

% Plot each group
groups = {'f';'m'};

for ng = 1:numel(groups)
    thisIds = eval(sprintf('ids_%s',groups{ng}));
    thisFits = eval(sprintf('fits_%s',groups{ng}));
    if ~contains(groups{ng}, 'm')
        thisC = plotParams.femaleC;
    else
        thisC = plotParams.maleC;
    end
    % Plot PDFs
    f = figure('Position',get(0,'ScreenSize'));
    for na = 1:numel(thisIds)
        % Plot histograms
        
        %%% Plot hisogram of non-laser trials
        subplot(size(thisFits,1)*2,round(numel(thisIds)/2)+1,na); box off; hold on
        data = Latency(strcmp(behaviorTable.aID,thisIds{na}));
        thishist = histcounts(data,bins,'normalization','pdf');
        sh=stairs(bins(1:end-1),thishist,'Color','none');
        % fill stairs
        X_plot = [sh.XData(1),repelem(sh.XData(2:end),2)];
        y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
        % Fill area
        fill([X_plot,fliplr(X_plot)],[y,zeros(size(y))],thisC,'LineStyle','none','FaceAlpha',.4)
        set(gca,'XLim',[bins(1),bins(end-1)],'YLim',[0 1])
        axis square
        
        %%% Plot fits
        
        y = custpdf(bins2,thisFits{na}(1), thisFits{na}(2), thisFits{na}(3));
        plot(bins2,y,':','Color',thisC,'LineWidth',1)
        axis square
        set(gca,'XLim',[bins(1),bins(end-1)]);%,'YLim',[0 1])
    end
end

%% Plot parameter means 
groups = {'f';'m'};

xaxis = repmat((1:2:numel(groups)*2)',1,2) + repmat([-.25,.25], numel(groups), 1);

f=figure('Position',[440 365 592 433]);

for ng = 1:numel(groups)
    
    thisIds = eval(sprintf('ids_%s',groups{ng}));
    thisFits = eval(sprintf('fits_%s',groups{ng}));
    if ~contains(groups{ng}, 'm') 
        thisC = plotParams.femaleC;
    else
        thisC = plotParams.maleC;
    end
        
    for nc = 1:size(thisFits,1) % for each condition, if relevant 
        clear X X_laser
        for na = 1:numel(thisIds)
            X(na,:) = thisFits{nc,na};
        end
        
        for np = 1:numParams
            p(np,nc) = subplot(size(thisFits,1),numParams,np+numParams*(nc-1)); hold on
            scatter(ones(size(X,1),1).*xaxis(ng,1),X(:,np),20,'MarkerFaceColor',thisC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
            plot([xaxis(ng,1)-.1 xaxis(ng,1)+.1] ,[mean(X(:,np)) mean(X(:,np))],'Color',[.3 .3 .3],'MarkerFaceColor',thisC,'MarkerEdgeColor','none');
            errorbar(xaxis(ng,1),mean(X(:,np)),std(X(:,np))./sqrt(size(X,1)),'LineStyle','none','Color',thisC,'LineWidth',1,'CapSize',0)
            ylabel(sprintf('\\%s',thisParams{np}));
            axis square
        end
        eval(sprintf('stats.ids_%s(:,nc) = thisIds;', groups{ng}));
        eval(sprintf('stats.param_%s{nc} = X;', groups{ng}));
    end
end

for nc = 1:size(thisFits,1)
    for np = 1:numParams
        set(p(np,nc),'XTick',1:2:numel(groups)*2,'XTickLabel',groups,'XTickLabelRotation',45,'XLim',[0 .5+numel(groups)*2]);
        if size(thisFits,1)>1
            set(p(np,nc),'YLim',[min([p(np,1).YLim(1) p(np,2).YLim(1)]) max([p(np,1).YLim(2) p(np,2).YLim(2)])]);
        end
    end
end

stats.fits_f = fits_f;
stats.fits_m = fits_m; 
  
%% Stats
groups = {'f';'m'};

vec = cell(1,numParams);
groupvec = cell(1,numParams);
ids_f = 1:numel(ids_f);
ids_m = 1+numel(ids_f):numel(ids_m)+numel(ids_f);
for ng = 1:numel(groups)
    thisIds = eval(sprintf('ids_%s',groups{ng}));
    thisFits = eval(sprintf('fits_%s',groups{ng}));
    clear X X_laser
    
    for nc = 1:size(thisFits,1)
        for na = 1:numel(thisIds)
            X(na,:,nc) = thisFits{nc,na};
        end
    end
    
    for np = 1:numParams
        for nc = 1:size(thisFits,1)
            vec{np} = cat(1,vec{np},X(:,np,nc));
            if size(thisFits,1)>1
                groupvec{np} = cat(1,groupvec{np},...
                    cat(2,ones(size(thisIds')).*~contains(groups{ng},'m'),... % sex F=1 M=0
                    thisIds',... % Animal ID
                    ones(size(thisIds)).*(nc-1))); % groupingVar 1
                
            else
                groupvec{np} = cat(1,groupvec{np},...
                    cat(2,ones(size(thisIds')).*~contains(groups{ng},'m'),... % sex F=1 M=0
                    thisIds')); %Animal ID 
            end
        end
    end
end
% 


if size(thisFits,1)==1
   for np = 1:numParams     
      [stats.ranksum.p(np),~,stats.ranksum.stats{np}] =  ranksum(vec{np}(groupvec{np}(:,1)==0),vec{np}(groupvec{np}(:,1)==1));
   end
end

