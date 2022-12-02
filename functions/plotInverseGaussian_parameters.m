function stats = plotInverseGaussian_parameters(fits,groups)

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

if nargin < 2
    groups = {'f';'m';'f_yfp';'m_yfp'};
end

numParams = numel(fits.thisParams); 

%% Plot parameter means 

xaxis = repmat((1:2:numel(groups)*2)',1,2) + repmat([-.25,.25], numel(groups), 1);

f=figure('Position',[440 654 364 134]);

for ng = 1:numel(groups)
    
    thisIds = eval(sprintf('fits.ids_%s',groups{ng}));
    thisFits = eval(sprintf('fits.fits_%s',groups{ng}));
    if ~contains(groups{ng}, 'm') 
        thisC = plotParams.femaleC;
    else
        thisC = plotParams.maleC;
    end
        
    for nc = 1:size(thisFits,1) % for each condition, if relevant 
        clear X X_laser
        for na = 1:numel(thisIds)
            X(na,:) = thisFits{nc,na}{1};
            try
            X_laser(na,:) = thisFits{nc,na}{2};
            catch
                keyboard
            end
        end
        
        for np = 1:numParams
            p(np,nc) = subplot(size(thisFits,1),numParams,np+numParams*(nc-1)); hold on
            
            scatter(ones(size(X,1),1).*xaxis(ng,1),X(:,np),15,'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
            scatter(ones(size(X_laser,1),1).*xaxis(ng,2),X_laser(:,np),15,'MarkerFaceColor',thisC,'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
            
            plot(xaxis(ng,1),nanmean(X(:,np)),'o','Color',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none');
            errorbar(xaxis(ng,1),nanmean(X(:,np)),nanstd(X(:,np))./sqrt(sum(~isnan(X_laser(:,np)))),'LineStyle','none','Color',[.3 .3 .3],'LineWidth',1)
            plot(xaxis(ng,2),nanmean(X_laser(:,np)),'o','Color',thisC,'MarkerFaceColor',thisC,'MarkerEdgeColor','none');
            errorbar(xaxis(ng,2),nanmean(X_laser(:,np)),nanstd(X_laser(:,np))./sqrt(sum(~isnan(X_laser(:,np)))),'LineStyle','none','Color',thisC,'LineWidth',1)
            ylabel(sprintf('\\%s',fits.thisParams{np}));
        end
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

  
%% Stats
groups = {'f';'m';'f_yfp';'m_yfp'};

vec = cell(1,numParams);
groupvec = cell(1,numParams);
vec2 = cell(1,numParams);
groupvec2 = cell(1,numParams);

for ng = 1:numel(groups)
    thisIds = eval(sprintf('fits.ids_%s',groups{ng}));
    thisFits = eval(sprintf('fits.fits_%s',groups{ng}));
    clear X X_laser
    
    for nc = 1:size(thisFits,1)
        for na = 1:numel(thisIds)
            X(na,:,nc) = thisFits{nc,na}{1};
            X_laser(na,:,nc) = thisFits{nc,na}{2};
        end
    end
    eval(sprintf('X_%s = X;', groups{ng}));
    eval(sprintf('X_laser_%s = X_laser;',groups{ng}))
    
    for np = 1:numParams
        for nc = 1:size(thisFits,1)
            % Laser vs. non-laser
            eval(sprintf('[stats.%s_p(np,nc),~,stats.%s_stats{np,nc}] = signrank(X(:,np,nc),X_laser(:,np,nc));',groups{ng},groups{ng}));
        end
    end
end


