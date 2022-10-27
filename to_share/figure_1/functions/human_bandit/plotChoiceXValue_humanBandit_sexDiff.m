function plotChoiceXValue_humanBandit_sexDiff(choice_f,choice_m,bins,valType)     

load(fullfile(whereAreWe('bucket'), 'Manuscript_figures','plotParams.mat'))


for nv = 1:numel(valType)
    f=figure('Units','inches','Position',[5,5,7,5]); hold on
    thisbins = eval(sprintf('bins.%s;',valType{nv}));
    thisbins = thisbins+mean(diff(thisbins))/2; % center bins
    thisplot = eval(sprintf('choice_f.probRight_%s;',valType{nv}));
    p=errorbar(thisbins,nanmean(thisplot),nanstd(thisplot)./sqrt(size(thisbins,1)),'Color',femaleC,'CapSize',0,'LineWidth',2);
    plot(repmat(thisbins,size(thisplot,1),1)',thisplot','Color',[femaleC .5],'LineWidth',.25);
    
    thisplot = eval(sprintf('choice_m.probRight_%s;',valType{nv}));
    p(2)=errorbar(thisbins,nanmean(thisplot),nanstd(thisplot)./sqrt(size(thisbins,1)),'Color',maleC,'CapSize',0,'LineWidth',2);
    plot(repmat(thisbins,size(thisplot,1),1)',thisplot','Color',[maleC .5],'LineWidth',.25);
    
    ylabel('P(choice=right)')
    
    switch valType{nv}
        case 'QTot'
            xlabel('Total value')
        case 'QChosenDiff'
            xlabel('Q_{chosen} - Q_{unchosen}')
        case 'QChosen'
            xlabel('Q_chosen')
        case 'QUnchosen'
            xlabel('Q_unchosen')
        case 'QDiff'
            xlabel('Q_{right} - Q_{left}')
    end
    
    set(gca,'XLim', [thisbins(1) thisbins(end)]-mean(diff(thisbins))/2,'FontSize',12)
    legend(p,{'Female','Male'})
    
    
end

