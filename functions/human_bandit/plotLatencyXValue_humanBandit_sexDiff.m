function stats=plotLatencyXValue_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag)

load(fullfile(whereAreWe('data'),'plotParams.mat'))


for nl = 1:numel(latencyType)
    for nv = 1:numel(valType)
        f=figure('Units','inches','Position',[5,5,7,5]); hold on
        thisbins = eval(sprintf('bins.%s;',valType{nv}));
        thisbins = 1:numel(thisbins); % center bins
        thisplot = eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv}));
        % plot(repmat(thisbins,size(thisplot,1),1)',thisplot','Color',[femaleC .25],'LineWidth',.25);
        p=errorbar(thisbins(1:end-1),nanmean(thisplot),nanstd(thisplot)./sqrt(size(thisplot,1)),'Color',femaleC,'CapSize',0,'LineWidth',2.5);
        
        thisplot = eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv}));
        % plot(repmat(thisbins,size(thisplot,1),1)',thisplot','Color',[maleC .25],'LineWidth',.25);
        p(2)=errorbar(thisbins(1:end-1),nanmean(thisplot),nanstd(thisplot)./sqrt(size(thisplot,1)),'Color',maleC,'CapSize',0,'LineWidth',2.5);
        
        
        if contains(latencyType{nl},'choice')
            if zscoreFlag
                ylabel('Choice latency (zscore)')
            else
                ylabel('Choice latency (ms)')
            end
        elseif contains(latencyType{nl},'trialStart')
            if zscoreFlag
                ylabel('Trial start latency (zscore)')
            else
                ylabel('Trial start latency (ms)')
            end
        end
        
        switch valType{nv}
            case 'QTot'
                xlabel('Total value quantiles')
            case 'QChosenDiff'
                xlabel('Q_{chosen} - Q_{unchosen} quantiles')
            case 'QChosen'
                xlabel('Q_chosen quantiles')
            case 'QUnchosen'
                xlabel('Q_unchosen quantiles')
            case 'QDiff'
                xlabel('Q_{right} - Q_{left} quantiles')
        end
        
        set(gca,'XLim', [thisbins(1)-.5 thisbins(end-1)+.5],'FontSize',12,'XTick',[thisbins(1):thisbins(end-1)])
        legend(p,{'Female','Male'})
        
        
        
        
        groupvec = [];
        vec = [];
        
        for nb = 1:numel(thisbins)-1
            vec=cat(1,vec,eval(sprintf('latency_m.%s_%s_quant(:,nb);',latencyType{nl},valType{nv})));
            groupvec = cat(1,groupvec,cat(2,ones(size(eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv})),1),1).*nb,ones(size(eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv})),1),1),latency_m.age',[1:numel(latency_m.age)]'));
            vec=cat(1,vec,eval(sprintf('latency_f.%s_%s_quant(:,nb);',latencyType{nl},valType{nv})));
            groupvec = cat(1,groupvec,cat(2,ones(size(eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv})),1),1).*nb,zeros(size(eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv})),1),1),latency_f.age',[1+numel(latency_m.age):numel(latency_m.age)+numel(latency_f.age)]'));
        end
        T = table(vec,'VariableNames',{'Latency'});
        T.value = categorical(groupvec(:,1));
        T.sex = categorical(groupvec(:,2));
        T.age = (groupvec(:,3));
        T.subject = categorical(groupvec(:,4));
        f = 'Latency ~ value*sex*age + (1|subject)';
        glme = fitlme(T,f,'DummyVarCoding','effects');
        stats{nv}.glme_quant = glme;
        stats{nv}.glme_quant_anova = dataset2table(anova(glme,'DFMethod','Satterthwaite'));
        for nb = 1:numel(thisbins)-1
            [stats{nv}.postHoc_rankSum_quant.p(nb),~,stats{nv}.postHoc_rankSum_quant.stats{nb}] = ranksum(vec(groupvec(:,2)==1&groupvec(:,1)==nb),vec(groupvec(:,2)==0&groupvec(:,1)==nb));
        end
        
        
        
    end
end