function stats = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins,behaviorTable,aids_f,aids_m,runReg)

%%% Figure 1h

% Plot latency x value x sex 

% latency_f and latency_m: average latencies by condition
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% latencyType: trialInit, leverPress
% bins: value bins
% behaviorTable: trial by trial data
% aids_f: list of female mice
% aids_m: list of male mice 

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
behaviorTable = behaviorTable(behaviorTable.laserSession==0,:); 
%% Plot trial initiation latency for males and females as a function of value (quantiles)
% Figure 1h
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
        eval(sprintf('stats.posthoc.%s_%s_quant = tempStats;', valType{nv}, latencyType{nl}))
    end
end

%% Mixed-effects model with weight, value, trial number and sex
if runReg
    for nv = 1:numel(valType)
        % select animals
        idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
        idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
        X = behaviorTable(cell2mat(idx),:);
        % select non-laser sessions
        X = X(X.laserSession==0,:);
        X.trial = nanzscore(X.trial);
        X.female = categorical(X.female);
        eval(sprintf('X.%s_quant = categorical(X.%s_quant);',valType{nv},valType{nv}));
        
        f = sprintf('trialInit_thresh ~ %s_quant*female*weight_zscore + %s_quant*female*trial + (1+%s_quant*trial + %s_quant*weight_zscore|aID)+(1+%s_quant*trial|aID:session)',valType{nv},valType{nv},valType{nv},valType{nv},valType{nv});
        
        mdl = fitlme(X,f,'DummyVarCoding','effects');
        stats.mdl{nv}.mdl = mdl;
        stats.mdl{nv}.mdl_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite'));
        stats.mdl{nv}.mdl_coeff = dataset2cell(mdl.Coefficients);
    end
end