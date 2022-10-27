function plotBilinearRegression(fit,r2,savename,savehere,fh)
% generate estimated data
X = [];
for ne = 1:numel(fit.eventNames)
    clear thisPred;
    clear coeff
    for ng = 1:numel(fit.gainNames{ne})+1
        coeff(:,ng) = fit.model_parameters.trial_coefficients{ne}(ng).*fit.gain_preds{ne,ng}(~isnan(fit.event_times{ne}));
    end
    coeff = sum(coeff,2);
    thisPred = fit.kernelPreds{ne};
    thisPred(thisPred==1) = coeff;
 
    for nb = 1:fit.nbasis(ne)
        thisbs = (full(fit.bs{ne}(:,nb)));
        temp = conv(thisPred,thisbs);
        X = cat(2,X,temp(1:numel(fit.kernelPreds{ne})));
    end
    
end
%X = nanzscore(X);
predY = X*cat(1,fit.model_parameters.kernel_profiles{:});

for np = 1:numel(fit.eventNames)
    %extract estimated data by event and gain predictor
    clear thistrials
    thisEvents = find(fit.kernelPreds{np}==1)+fit.eventDur(np,1)*fit.frameRate;
    for ne = 1:numel(thisEvents)
        try
            thistrials(ne,:) = predY(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
        catch
            thistrials(ne,:) = predY(1,predY(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:end),zeros(thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate-numel(predY),1));
        end
    end
    
    predY_event{np} = nanmean(thistrials);
    semPredY_event{np} = nanstd(thistrials)./sqrt(size(thistrials,1));
    for ng = 1:numel(fit.gainNames{np})
        [~,~,bins] = histcounts(fit.gain_preds{np,ng},fit.gainBins{ng});
        bins = bins(~isnan(fit.event_times{np}));
        for nb = 1:numel(fit.gainBins{ng})
            predFluor_bins{np,nb,ng} = nanmean(thistrials(bins==nb,:),1);
            %predFluor_value{np,nb,ng} = nanmean(thistrials(fit.gain_preds{np}==fit.gainBins{ng}(nb),:),1);
            semPredFluor_bins{np,nb,ng} = nanstd(thistrials(bins==nb,:),[],1)./sqrt(sum(bins==nb));
           % semPredFluor_value{np,nb,ng} = nanstd(thistrials(fit.gain_preds{np}==fit.gainBins{ng}(nb),:),[],1)./sqrt(sum(fit.gainPreds{np}==fit.gainBins{ng}(nb)));
        end
    end
end


% find kernels
for np = 1:numel(fit.eventNames)
    clear coeff
    for ng = 1:numel(fit.gainNames{np})
        coeff(:,ng) = fit.model_parameters.trial_coefficients{np}(ng).*fit.gain_preds{np,ng}(~isnan(fit.event_times{np}));
    end
    %coeff = nanmean(sum(coeff,2));
    coeff = 1;
    thisbs = (fit.bs{np});
    thisKernel{np} = (thisbs*fit.model_parameters.kernel_profiles{np}).*coeff;
    clear coeff
    
end
ylims = [min(min(cell2mat(thisKernel'))) max(max(cell2mat(thisKernel')))];


for np = 1:numel(fit.eventNames)
    nplots = (2*numel(fit.gainNames{np})+2)/2;
    figure(fh(1)); set(gcf,'Visible','off') 
    xaxis  = linspace(-fit.eventDur(np,1), fit.eventDur(np,2), size(fit.meanFluor_event{np},2));
    % plot trial-averaged predicted data
    subplot(nplots,2,1)
    hold off
    pp(1)=shadedErrorBar(xaxis,fit.meanFluor_event{np}, fit.semFluor_event{np});
    hold on
    pp(2)=shadedErrorBar(xaxis,predY_event{np},semPredY_event{np},'lineProps',{'--'});
    
    legend([pp(1).mainLine pp(2).mainLine], {'real';'predicted'},'Box','off')
    title(sprintf('%s \n r^2: %s',fit.eventNames{np},num2str(r2)));
    
    % Plot kernel
    subplot(nplots,2,2)  
    hold off
    plot(xaxis,thisKernel{np})
    hold on
    box off
    set(gca,'YLim',ylims)
    % plot estimated and real data by gain for each gain
    for ng = 1:numel(fit.gainNames{np})
        % Plot average activity by gain
        subplot(nplots,2,2*ng+1); hold off
        for nb = 1:numel(fit.gainBins{np,ng})-1
            temp = shadedErrorBar(xaxis,fit.meanFluor_bins{np,nb,ng},fit.semFluor_bins{np,nb,ng},'lineProps',{'LineWidth',1});
            hold on
            p(nb) = temp.mainLine;
            %shadedErrorBar(xaxis,predFluor_bins{np,nb,ng},semPredFluor_bins{np,nb,ng},'lineProps',{'--','LineWidth',1,'Color',p(nb).Color});
            title(sprintf('%s: %s',fit.gainNames{np}{ng},num2str(fit.model_parameters.trial_coefficients{np}(ng))))
        end
        legend(p,{num2str(fit.gainBins{np,ng}(1:end-1)')},'Box','off');
        box off
    end
    
    xlabel('Time from event (sec)')
    print(fullfile(savehere, sprintf('%s_%s',savename,fit.eventNames{np})), '-dpdf', '-bestfit')
end

figure(fh(2)); set(gcf,'Visible','off','WindowStyle','docked') 
nplots = numel(fit.eventNames);
for np = 1:nplots
    subplot(2,nplots,np)
    xaxis  = linspace(-fit.eventDur(np,1), fit.eventDur(np,2), size(fit.meanFluor_event{np},2));
    % plot trial-averaged predicted data
    hold off
    pp(1)=shadedErrorBar(xaxis,fit.meanFluor_event{np}, fit.semFluor_event{np});
  %  pp(2)=shadedErrorBar(xaxis,predY_event{np},semPredY_event{np},'lineProps',{'--'});
  hold on
    if np==1
    legend([pp(1).mainLine pp(2).mainLine], {'real';'predicted'},'Box','off')
    end
    title(sprintf('%s \n r^2: %s',fit.eventNames{np},num2str(r2)));
    
    % Plot kernel
    subplot(2,nplots,np+nplots)
    plot(xaxis,thisKernel{np})
    box off
    set(gca,'YLim',ylims)
end
print(fullfile(savehere, sprintf('%s_kernels',savename)), '-dpdf', '-bestfit')
