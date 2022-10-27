function fit = extractEvents_bilinear(fit,y,ptiles,stay)
if nargin<3
ptiles = [0 25 50 75 100];
end
for ne = 1:numel(fit.eventNames)
    for ng = 1:numel(fit.gainNames{ne})
        fit.gainBins{ne,ng} = prctile(fit.gain_preds{ne,ng},ptiles);
        %fit.gainBins{ne,ng} = min(fit.gain_preds{ne,ng}):.4:max(fit.gain_preds{ne,ng});
    end
end

fit.meanFluor_bins = [];
fit.semFluor_bins = [];
fit.meanFluor_value = [];
fit.semFluor_value = [];
for np = 1:numel(fit.event_times)
    try
    y_resid = extractResidual(np,fit,y);
    catch
        y_resid = nan(size(y));
    end
    clear thistrials thistrials_resid
    thisEvents = fit.event_times{np}';
    if size(thisEvents,2)>size(thisEvents,1)
        thisEvents = thisEvents';
    end
    for ne = 1:numel(thisEvents)
        if isnan(thisEvents(ne))
            thistrials(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
            thistrials_resid(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
        else
            try
                thistrials(ne,:) = y(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
                thistrials_resid(ne,:) = y_resid(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
            catch
                thistrials(ne,:) = cat(1,y(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:end),zeros(thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate-numel(y),1));
                thistrials_resid(ne,:) = cat(1,y_resid(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:end),zeros(thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate-numel(y),1));
            end
        end
    end
    fit.meanFluor_event{np} = nanmean(thistrials);
    fit.semFluor_event{np} = nanstd(thistrials)./sqrt(sum(~isnan(thisEvents)));
    fit.meanResid_event{np} = nanmean(thistrials_resid); 
    fit.semResid_event{np} = nanstd(thistrials_resid)./sqrt(sum(~isnan(thisEvents)));

    for ng = 1:numel(fit.gainNames{np})
        [~,~,bins] = histcounts(fit.gain_preds{np,ng},fit.gainBins{np,ng});
        %[~,~,bins] = histcounts(fit.gain_preds{np,ng},linspace(min(fit.gain_preds{np,ng}),max(fit.gain_preds{np,ng}),numel(ptiles)));
        for nb = 1:numel(fit.gainBins{np,ng})
           
            fit.meanFluor_bins{np,nb,ng} = nanmean(thistrials(bins==nb,:),1);
%            if isnan(fit.meanFluor_bins{np,nb,ng}) & nb~=numel(fit.gainBins{np,ng})
%                keyboard
%            end
            fit.meanFluor_value{np,nb,ng} = nanmean(thistrials(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),1);
            
            fit.meanResid_bins{np,nb,ng} = nanmean(thistrials_resid(bins==nb,:),1);
            fit.meanResid_value{np,nb,ng} = nanmean(thistrials_resid(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),1);
            
            fit.semFluor_bins{np,nb,ng} = nanstd(thistrials(bins==nb,:),[],1)./sqrt(sum(bins==nb&~isnan(thisEvents)));
            fit.semFluor_value{np,nb,ng} = nanstd(thistrials(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),[],1)./sqrt(sum(~isnan(thisEvents) & fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb)));
            fit.semResid_bins{np,nb,ng} = nanstd(thistrials_resid(bins==nb,:),[],1)./sqrt(sum(bins==nb&~isnan(thisEvents)));
            fit.semResid_value{np,nb,ng} = nanstd(thistrials_resid(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),[],1)./sqrt(sum(~isnan(thisEvents) & fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb)));

        end
    end
    
    
    if nargin == 4
        for ng = 1:numel(fit.gainNames{np})
            [~,~,bins] = histcounts(fit.gain_preds{np,ng},fit.gainBins{np,ng});
            for nb = 1:numel(fit.gainBins{np,ng})
                fit.meanFluor_stay_bins{np,nb,ng} = nanmean(thistrials(bins==nb&stay==1,:),1);
                fit.semFluor_stay_bins{np,nb,ng} = nanstd(thistrials(bins==nb&stay==1,:),[],1)./sqrt(sum(bins==nb&stay==1&~isnan(thisEvents)));
                fit.meanFluor_switch_bins{np,nb,ng} = nanmean(thistrials(bins==nb&stay==0,:),1);
                fit.semFluor_switch_bins{np,nb,ng} = nanstd(thistrials(bins==nb&stay==0,:),[],1)./sqrt(sum(bins==nb&stay==1&~isnan(thisEvents)));
            end
        end 
    end
    
    
    
end
end

function y_resid = extractResidual(whichEvent,fit,y)
X = [];

for ne = setdiff(1:numel(fit.eventNames),whichEvent)
    clear coeff
    for ng = 1:numel(fit.gainNames{ne})+1
        coeff(:,ng) = fit.model_parameters.trial_coefficients{ne}(ng).*fit.gain_preds{ne,ng}(~isnan(fit.event_times{ne}));
    end
    coeff = sum(coeff,2);
    
    thisPred = fit.kernelPreds{ne};
    thisPred(thisPred==1) = coeff;
    
    for nb = 1:fit.nbasis(ne)
        temp = conv(thisPred,full(fit.bs{ne}(:,nb)));
        X = cat(2,X,temp(1:numel(y)));
    end
end
y_pred = X*cat(1,fit.model_parameters.kernel_profiles{setdiff(1:numel(fit.eventNames),whichEvent)}); % predicted values
y_resid = y-y_pred;
end

    
    
