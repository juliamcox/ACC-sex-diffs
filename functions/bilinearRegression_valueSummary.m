function bilinearRegression_valueSummary(regVer,dataset,ptiles,toPlot,rawFlag,frameRate,qFile)
if nargin < 5
    rawFlag = 1;
    frameRate = 20;
    zscoreGain = 0;
    qFile = 'qLearn_session_all';
end
extractEvents = 1;
a=0.01;

[events, kernelType, gains, eventDur,nbasis] = getEvents_bilinear(regVer);

fbasename = fullfile(whereAreWe('bilinear'),dataset,regVer);

try load(fullfile(fbasename,'summary.mat'))
catch
    bilinearRegression_summary(regVer,dataset,rawFlag)
    load(fullfile(fbasename,'summary.mat'))
end

% find significant neurons
pmat_f = pvals<a;
pmat_m = pvals_m<a;

eventIdx = 1;
for ne = 1:numel(events)
    for ng = 1:numel(gains{ne})+1
        sigNeurons_f{ne,ng} = find(pmat_f(:,eventIdx));
        sigNeurons_m{ne,ng} = find(pmat_m(:,eventIdx));
        eventIdx = eventIdx+1;
    end
end

try
load(fullfile(whereAreWe('bilinear'),sprintf('%s_activityXvalue.mat',regVer))); 

catch
%% Plot trial-averaged activity and predicted activity

% Load each neuron's data and extract predicted and real data for plotting

for ne = 1:numel(events)
    %   keyboard
    for ng = 1:numel(gains{ne})
        predPlot_neg_f{ne,ng} = [];
        predPlot_neg_m{ne,ng} = [];
        predPlot_pos_f{ne,ng} = [];
        predPlot_pos_m{ne,ng} = [];
        realPlot_neg_f{ne,ng} = [];
        realPlot_neg_m{ne,ng} = [];
        realPlot_pos_f{ne,ng} = [];
        realPlot_pos_m{ne,ng} = [];
        realPlot_f{ne,ng} = [];
        realPlot_m{ne,ng} = [];
        predPlot_f{ne,ng} = [];
        predPlot_m{ne,ng} = [];
        for nn = 1:numel(sigNeurons_f{ne,ng})
            load(fullfile(whereAreWe('bilinear'),dataset,sprintf('rawFlag%d_frameRate%d_%s',rawFlag,frameRate,qFile),thisFnames{sigNeurons_f{ne,ng}(nn)}));
            load(fullfile(fbasename,thisFnames{sigNeurons_f{ne,ng}(nn)}),'fit_real');
            
            if extractEvents
                fit_real = extractEvents_bilinear(fit_real,nanzscore(y),ptiles);
                [predY_value, ~] = calcPredData(fit_real,y);
            else
                try
                    [predY_value, ~] = calcPredData(fit_real,y);
                catch ME
                    keyboard
                end
            end
            
            thisKernel = fit_real.bs{ne}*fit_real.model_parameters.kernel_profiles{ne};
            if trapz(thisKernel) > 0 && trapz(fit_real.meanFluor_event{ne}) > 0  % if kernel and real data are the same sign
                if fit_real.model_parameters.trial_coefficients{ne}(ng)<0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_neg_f{ne,ng} = cat(1,realPlot_neg_f{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_neg_f{ne,ng} = cat(1,predPlot_neg_f{ne,ng}, temp);
                elseif fit_real.model_parameters.trial_coefficients{ne}(ng)>0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    try
                        realPlot_pos_f{ne,ng} = cat(1,realPlot_pos_f{ne,ng}, temp);
                    catch
                        keyboard
                    end
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_pos_f{ne,ng} = cat(1,predPlot_pos_f{ne,ng}, temp);
                else
                    keyboard
                end
            elseif trapz(thisKernel) < 0 && trapz(fit_real.meanFluor_event{ne}) > 0 
                if fit_real.model_parameters.trial_coefficients{ne}(ng)>0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_neg_f{ne,ng} = cat(1,realPlot_neg_f{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_neg_f{ne,ng} = cat(1,predPlot_neg_f{ne,ng}, temp);
                elseif fit_real.model_parameters.trial_coefficients{ne}(ng)<0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_pos_f{ne,ng} = cat(1,realPlot_pos_f{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_pos_f{ne,ng} = cat(1,predPlot_pos_f{ne,ng}, temp);
                else
                    keyboard
                end
            else
                clear temp
                temp = fit_real.meanFluor_bins(ne,:,ng);
                temp = cell2mat(temp(1:end-1)');
                temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                temp = permute(temp,[2 1 3]);
                
                realPlot_neg_f{ne,ng} = cat(1,realPlot_neg_f{ne,ng}, temp);
                clear temp
                temp = predY_value(ne,:,ng);
                temp = cell2mat(temp(1:end-1)');
                temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                temp = permute(temp,[2 1 3]);
                predPlot_neg_f{ne,ng} = cat(1,predPlot_neg_f{ne,ng}, temp);
            end
            clear temp
            temp = fit_real.meanFluor_bins(ne,:,ng);
            temp = cell2mat(temp(1:end-1)');
            temp = (reshape(temp',size(temp,2),1,size(temp,1)));
            temp = permute(temp,[2 1 3]);
            realPlot_f{ne,ng} = cat(1,realPlot_f{ne,ng},temp);
            clear temp
            temp = predY_value(ne,:,ng);
            temp = cell2mat(temp(1:end-1)');
            temp = (reshape(temp',size(temp,2),1,size(temp,1)));
            temp = permute(temp,[2 1 3]);
            predPlot_f{ne,ng} = cat(1,predPlot_f{ne,ng},temp);
        end
        
        for nn = 1:numel(sigNeurons_m{ne,ng})
            
            load(fullfile(fullfile(whereAreWe('bilinear'),dataset,sprintf('rawFlag%d_frameRate%d_%s',rawFlag,frameRate,qFile),thisFnames_m{sigNeurons_m{ne,ng}(nn)})));
            load(fullfile(fbasename,thisFnames_m{sigNeurons_m{ne,ng}(nn)}),'fit_real');
            if extractEvents
                fit_real = extractEvents_bilinear(fit_real,nanzscore(y),ptiles);
            end
            try
                [predY_value, ~] = calcPredData(fit_real,y);
            catch ME
                keyboard
            end
            thisKernel = fit_real.bs{ne}*fit_real.model_parameters.kernel_profiles{ne};
            if corr(thisKernel,fit_real.meanFluor_event{ne}')>0
                if fit_real.model_parameters.trial_coefficients{ne}(ng)<0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_neg_m{ne,ng} = cat(1,realPlot_neg_m{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_neg_m{ne,ng} = cat(1,predPlot_neg_m{ne,ng}, temp);
                elseif fit_real.model_parameters.trial_coefficients{ne}(ng)>0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_pos_m{ne,ng} = cat(1,realPlot_pos_m{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_pos_m{ne,ng} = cat(1,predPlot_pos_m{ne,ng}, temp);
                else
                    keyboard
                end
            elseif corr(thisKernel,fit_real.meanFluor_event{ne}')<0 % if the kernel is flipped, flip the sign of the gain coefficient
                if fit_real.model_parameters.trial_coefficients{ne}(ng)>0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_neg_m{ne,ng} = cat(1,realPlot_neg_m{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_neg_m{ne,ng} = cat(1,predPlot_neg_m{ne,ng}, temp);
                elseif fit_real.model_parameters.trial_coefficients{ne}(ng)<0
                    clear temp
                    temp = fit_real.meanFluor_bins(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    realPlot_pos_m{ne,ng} = cat(1,realPlot_pos_m{ne,ng}, temp);
                    clear temp
                    temp = predY_value(ne,:,ng);
                    temp = cell2mat(temp(1:end-1)');
                    temp = (reshape(temp',size(temp,2),1,size(temp,1)));
                    temp = permute(temp,[2 1 3]);
                    predPlot_pos_m{ne,ng} = cat(1,predPlot_pos_m{ne,ng}, temp);
                else
                    keyboard
                end
            else
                keyboard
            end
            temp = fit_real.meanFluor_bins(ne,:,ng);
            temp = cell2mat(temp(1:end-1)');
            temp = (reshape(temp',size(temp,2),1,size(temp,1)));
            temp = permute(temp,[2 1 3]);
            realPlot_f{ne,ng} = cat(1,realPlot_f{ne,ng},temp);
            clear temp
            temp = predY_value(ne,:,ng);
            temp = cell2mat(temp(1:end-1)');
            temp = (reshape(temp',size(temp,2),1,size(temp,1)));
            temp = permute(temp,[2 1 3]);
            predPlot_f{ne,ng} = cat(1,predPlot_f{ne,ng},temp);
        end
    end
end
save(fullfile(whereAreWe('bilinear'),sprintf('%s_activityXvalue.mat',regVer)),'realPlot_neg_f','realPlot_neg_m','realPlot_pos_f','realPlot_pos_m');
end
for ne = toPlot
    %   keyboard
    for ng = 1:numel(gains{ne})
        %% Plot average activity of positive and negatively encoding neurons by bin
        plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
        [mmap, fmap]   = maleFemaleColormap(plotParams.maleC,plotParams.femaleC,size(realPlot_neg_f{ne,ng},3)+1);
        mmap = flipud(mmap);
        fmap = flipud(fmap);
        
        pos_m = nanmean(realPlot_pos_m{ne,ng},2);
        neg_m = nanmean(realPlot_neg_m{ne,ng},2);
        
        thisPlot_m = cat(1,reshape(neg_m,size(neg_m,1),size(neg_m,3)), reshape(pos_m,size(pos_m,1),size(pos_m,3)));
        
        pos_f = nanmean(realPlot_pos_f{ne,ng},2);
        neg_f = nanmean(realPlot_neg_f{ne,ng},2);
        
        if isempty(pos_f)
            thisPlot_f = reshape(neg_f,size(neg_f,1),size(neg_f,3));
        else
            thisPlot_f = cat(1,reshape(neg_f,size(neg_f,1),size(neg_f,3)), reshape(pos_f,size(pos_f,1),size(pos_f,3)));
        end
        
        
        figure();
        subplot(1,3,1);hold on
        plot(repmat(1:size(thisPlot_m,2),size(thisPlot_f,1),1)', thisPlot_f','Color',[plotParams.femaleC .7])
        set(gca,'YLim',[-.4 .8],'XTick',1:size(thisPlot_m,2),'XLim',[0.5 size(thisPlot_m,2)+.5])
        ylabel('Mean fluorescence')
        xlabel('value')
        title(gains{ne}{ng})
        subplot(1,3,2);hold on
        plot(repmat(1:size(thisPlot_m,2),size(thisPlot_m,1),1)', thisPlot_m','Color',[plotParams.maleC .7])
        set(gca,'YLim',[-.4 .8],'XTick',1:size(thisPlot_m,2),'XLim',[0.5 size(thisPlot_m,2)+.5])
        title(sprintf('n_f = %d, n_m = %d',size(thisPlot_f,1),size(thisPlot_m,1)))
        ylabel('Mean fluorescence')
        xlabel('value')
        
        subplot(1,3,3);hold on
        plot(nanmean(thisPlot_f,1),'Color',plotParams.femaleC, 'LineWidth',1.5)
        errorbar(1:size(thisPlot_f,2), nanmean(thisPlot_f,1), nanstd(thisPlot_f,[],1)./sqrt(size(thisPlot_f,1)),'LineStyle','none','LineWidth',1.5,'CapSize',0,'Color',plotParams.femaleC)
        plot(nanmean(thisPlot_m,1),'Color',plotParams.maleC, 'LineWidth',1.5)
        errorbar(1:size(thisPlot_f,2), nanmean(thisPlot_m,1), nanstd(thisPlot_m,[],1)./sqrt(size(thisPlot_m,1)),'LineStyle','none','LineWidth',1.5,'CapSize',0,'Color',plotParams.maleC)
        set(gca,'YLim',[-.15 .2],'XTick',1:size(thisPlot_m,2),'XLim',[0.5 size(thisPlot_m,2)+.5])
        ylabel('Mean fluorescence')
        xlabel('value')
        
        
        
    end
    
end
end



function [predFluor_bins, predFluor_value,predY_event] = calcPredData(fit,y)
y = nanzscore(y);
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
        try
        X = cat(2,X,temp(1:numel(fit.kernelPreds{ne})));
        catch
            keyboard
        end
    end
    
end
predY = X*cat(1,fit.model_parameters.kernel_profiles{:});

for np = 1:numel(fit.eventNames)
    %extract estimated data by event and gain predictor
    clear thistrials
    thisEvents = fit.event_times{np}';
    for ne = 1:numel(thisEvents)
        if isnan(thisEvents(ne))
            thistrials(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
            thistrials2(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
        else
            thistrials(ne,:) = predY(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
            thistrials2(ne,:) = y(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
        end
    end

    predY_event{np} = nanmean(thistrials);
    semPredY_event{np} = nanstd(thistrials)./sqrt(sum(~isnan(thisEvents)));
    for ng = 1:numel(fit.gainNames{np})
        [~,~,bins] = histcounts(fit.gain_preds{np,ng},fit.gainBins{np,ng});
        for nb = 1:numel(fit.gainBins{np,ng})
            predFluor_bins{np,nb,ng} = nanmean(thistrials(bins==nb,:),1);
            predFluor_value{np,nb,ng} = nanmean(thistrials(fit.gain_preds{np,ng}==fit.gainBins{ng}(nb),:),1);
            semPredFluor_bins{np,nb,ng} = nanstd(thistrials(bins==nb,:),[],1)./sqrt(sum(bins==nb));
            semPredFluor_value{np,nb,ng} = nanstd(thistrials(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),[],1)./sqrt(sum(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb)));

        end
    end
end


end