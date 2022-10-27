function bilinearRegression_valueSummary2(regVer,dataset,ptiles,toPlot,rawFlag,frameRate,zscoreGain,qFile)
if nargin < 5
    rawFlag = 1;
    frameRate = 20;
    zscoreGain = 0; 
    qFile = 'qLearn_session_all';
end
extractEvents = 0;
a=0.01;

[events, kernelType, gains, eventDur,nbasis] = getEvents_bilinear(regVer);
savehere = fullfile(whereAreWe('figurecode'),'processed_data','fig6'); 
savename = 'highLowValue';

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
            try
                load(fullfile(whereAreWe('bilinear'),dataset,regVer,thisFnames{sigNeurons_f{ne,ng}(nn)}),'fit_real');
            catch
                clear fit_real
                try
                load(fullfile(whereAreWe('bilinear'),dataset,regVer,thisFnames{sigNeurons_f{ne,ng}(nn)}),'fit_real');
                catch
                    keyboard
                end
            end
            if extractEvents
                load(fullfile(whereAreWe('imaging'),'bilinearRegression',dataset,sprintf('rawFlag%d_frameRate%d_%s',rawFlag,frameRate,qFile),thisFnames{sigNeurons_f{ne,ng}(nn)}));
                fit_real = extractEvents_bilinear(fit_real,nanzscore(y),ptiles);
                [predY_value, ~] = calcPredData(fit_real,y);
            end
%            [predY_value, ~] = calcPredData(fit_real);
            thisKernel = fit_real.bs{ne}*fit_real.model_parameters.kernel_profiles{ne};
            if trapz(thisKernel) > 0 && trapz(fit_real.meanFluor_event{ne}) > 0  % if kernel and real data are the same sign
                %corr(thisKernel,fit_real.meanFluor_event{ne}')>0
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
            elseif trapz(thisKernel) < 0 && trapz(fit_real.meanFluor_event{ne}) > 0 %|| trapz(thisKernel) < 0 && trapz(fit_real.meanFluor_events{ne}) < 0
                %corr(thisKernel,fit_real.meanFluor_event{ne}')<0 % if the kernel is flipped, flip the sign of the gain coefficient 
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
            try
            load(fullfile(whereAreWe('imaging'),'bilinearRegression',dataset,regVer,thisFnames_m{sigNeurons_m{ne,ng}(nn)}),'fit_real');
            catch
                clear fit_real
                try
            load(fullfile(whereAreWe('imaging'),'bilinearRegression',dataset,regVer,thisFnames_m{sigNeurons_m{ne,ng}(nn)}),'fit_real');
                catch
                    keyboard
                end
            end
            if extractEvents
                load(fullfile(whereAreWe('imaging'),'bilinearRegression',dataset,sprintf('rawFlag%d_frameRate%d_%s',rawFlag,frameRate,qFile),thisFnames_m{sigNeurons_m{ne,ng}(nn)}));
                fit_real = extractEvents_bilinear(fit_real,nanzscore(y),ptiles);
            end
            
            [predY_value, ~] = calcPredData(fit_real,y);
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

 for ne = toPlot
 %   keyboard
    for ng = 1:numel(gains{ne})
        %% Plot average activity of positive and negatively encoding neurons by bin
        plotParams = load(fullfile(whereAreWe('bucket'),'Manuscript_figures','plotParams.mat'));
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
        subplot(1,2,1);hold on
        plot(repmat(1:size(thisPlot_m,2),size(thisPlot_f,1),1)', thisPlot_f','Color',[plotParams.femaleC .7])
        plot(nanmean(thisPlot_f,1),'Color',plotParams.femaleC, 'LineWidth',2)
        set(gca,'YLim',[-.4 1.2],'XTick',1:size(thisPlot_m,2))
        subplot(1,2,2);hold on
        plot(repmat(1:size(thisPlot_m,2),size(thisPlot_m,1),1)', thisPlot_m','Color',[plotParams.maleC .7])
        plot(nanmean(thisPlot_m,1),'Color',plotParams.maleC, 'LineWidth',2)
        set(gca,'YLim',[-.4 1.2],'XTick',1:size(thisPlot_m,2))
        
        figure();
        subplot(1,2,1);hold on
        plot(nanmean(thisPlot_f,1),'Color',plotParams.femaleC, 'LineWidth',1.5)
        errorbar(1:size(thisPlot_f,2), nanmean(thisPlot_f,1), nanstd(thisPlot_f,[],1)./sqrt(size(thisPlot_f,1)),'LineStyle','none','LineWidth',1.5,'CapSize',0,'Color',plotParams.femaleC)
        plot(nanmean(thisPlot_m,1),'Color',plotParams.maleC, 'LineWidth',2)
        errorbar(1:size(thisPlot_f,2), nanmean(thisPlot_m,1), nanstd(thisPlot_m,[],1)./sqrt(size(thisPlot_m,1)),'LineStyle','none','LineWidth',1.5,'CapSize',0,'Color',plotParams.maleC)
        set(gca,'YLim',[-.15 .2],'XTick',1:size(thisPlot_m,2))
        
        negIdx = find(thisPlot_f(:,1)>thisPlot_f(:,end));
        posIdx = find(thisPlot_f(:,1)<thisPlot_f(:,end));
        negIdx_m = find(thisPlot_m(:,1)>thisPlot_m(:,end));
        posIdx_m = find(thisPlot_m(:,1)<thisPlot_m(:,end));
        
        
        
        figure()
        subplot(1,2,1); hold on
        errorbar(1:size(thisPlot_f(negIdx,:),2),nanmean(thisPlot_f(negIdx,:),1),std(thisPlot_f(negIdx,:))./sqrt(numel(negIdx)),'Color',plotParams.femaleC,'CapSize',0)
        errorbar(1:size(thisPlot_m(negIdx_m,:),2),nanmean(thisPlot_m(negIdx_m,:),1),std(thisPlot_m(negIdx_m,:))./sqrt(numel(negIdx_m)),'Color',plotParams.maleC,'CapSize',0)
        set(gca,'XTick',1:size(thisPlot_m,2),'XLim',[.5 size(thisPlot_m,2)+.5],'YLim',[-.4,1])
        xlabel('quantile')
        ylabel('nanmean fluorescence (z-score)')
        title(sprintf('n_f = %d, n_m = %d',numel(negIdx),numel(negIdx_m)))
        subplot(1,2,2); hold on
        errorbar(1:size(thisPlot_f(posIdx,:),2),nanmean(thisPlot_f(posIdx,:),1),std(thisPlot_f(posIdx,:),[],1)./sqrt(numel(posIdx)),'Color',plotParams.femaleC,'CapSize',0)
        errorbar(1:size(thisPlot_m(posIdx_m,:),2),nanmean(thisPlot_m(posIdx_m,:),1),std(thisPlot_m(posIdx_m,:))./sqrt(numel(posIdx_m)),'Color',plotParams.maleC,'CapSize',0)
        set(gca,'XTick',1:size(thisPlot_m,2),'XLim',[.5 size(thisPlot_m,2)+.5],'YLim',[-.4,1])
        xlabel('quantile')
        ylabel('nanmean fluorescence (z-score)')
        title(sprintf('n_f = %d, n_m = %d',numel(posIdx),numel(posIdx_m)))
        
        
        figure();
        subplot(2,2,1);hold on
        plot(repmat(1:size(thisPlot_m,2),numel(posIdx),1)', thisPlot_f(posIdx,:)','Color',[plotParams.femaleC .75],'LineWidth',.75)
        set(gca,'YLim',[-.5 2.5],'XTick',1:size(thisPlot_m,2))
        subplot(2,2,2);hold on
        plot(repmat(1:size(thisPlot_m,2),numel(posIdx_m),1)', thisPlot_m(posIdx_m,:)','Color',[plotParams.maleC .75],'LineWidth',.75)
        set(gca,'YLim',[-.5 2.5],'XTick',1:size(thisPlot_m,2))
        subplot(2,2,1);hold on
        plot(repmat(1:size(thisPlot_m,2),numel(negIdx),1)', thisPlot_f(negIdx,:)','--','Color',[plotParams.femaleC])
        set(gca,'YLim',[-.5 1],'XTick',1:size(thisPlot_m,2))
        set(gca,'XLim',[.75 3.25])
        subplot(2,2,2);hold on
        plot(repmat(1:size(thisPlot_m,2),numel(negIdx_m),1)', thisPlot_m(negIdx_m,:)','--','Color',[plotParams.maleC])
        set(gca,'YLim',[-.5 1],'XTick',1:size(thisPlot_m,2))
        set(gca,'XLim',[.75 3.25])
        
         %         figure()
%         subplot(1,2,1); hold on
%         errorbar(1:size(thisPlot_f,2),nanmean(thisPlot_f,1),std(thisPlot_f)./sqrt(size(thisPlot_f,1)),'Color',plotParams.femaleC,'CapSize',0)
%         errorbar(1:size(thisPlot_m,2),nanmean(thisPlot_m,1),std(thisPlot_m)./sqrt(size(thisPlot_m,1)),'Color',plotParams.maleC,'CapSize',0)
%         set(gca,'XTick',1:size(thisPlot_m,2),'XLim',[.5 size(thisPlot_m,2)+.5])
%         xlabel('quantile')
%         ylabel('nanmean fluorescence (z-score)')
        
%         realPlot_f = cat(1,realPlot_neg_f{ne,ng},realPlot_pos_f{ne,ng});
%         realPlot_neg_f{ne,ng} = realPlot_f(negIdx,:,:);
%         realPlot_pos_f{ne,ng} = realPlot_f(posIdx,:,:);
%         
%         realPlot_m = cat(1,realPlot_neg_m{ne,ng},realPlot_pos_m{ne,ng});
%         realPlot_neg_m{ne,ng} = realPlot_m(negIdx,:,:);
%         realPlot_pos_m{ne,ng} = realPlot_m(posIdx,:,:);
% 
%         xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_f{ne,ng},2));
%         figure(),
%         for np = 1:size(realPlot_neg_f{ne,ng},3)
%              subplot(2,2,1); hold on
%             shadedErrorBar(xaxis, nanmean(realPlot_neg_f{ne,ng}(:,:,np),1),nanstd(realPlot_neg_f{ne,ng}(:,:,np),[],1)./sqrt(size(realPlot_neg_f{ne,ng},1)),'lineprops',{'Color',fmap(np,:)})
%              subplot(2,2,2); hold on
%             shadedErrorBar(xaxis, nanmean(realPlot_neg_m{ne,ng}(:,:,np),1),nanstd(realPlot_neg_m{ne,ng}(:,:,np),[],1)./sqrt(size(realPlot_neg_m{ne,ng},1)),'lineprops',{'Color',mmap(np,:)})  
%             title('negative gain')
%         end
%         xlabel('Time from outcome (s)')
%         ylabel('nanmean z-scored fluorescence')
%         if size(realPlot_pos_f{ne,ng},1)~=0
%         for np = 1:size(realPlot_pos_f{ne,ng},3)
%              subplot(2,2,3); hold on
%             shadedErrorBar(xaxis, nanmean(realPlot_pos_f{ne,ng}(:,:,np),1),nanstd(realPlot_pos_f{ne,ng}(:,:,np),[],1)./sqrt(size(realPlot_pos_f{ne,ng},1)),'lineprops',{'Color',fmap(np,:)})
%              subplot(2,2,4); hold on
%             shadedErrorBar(xaxis, nanmean(realPlot_pos_m{ne,ng}(:,:,np),1),nanstd(realPlot_pos_m{ne,ng}(:,:,np),[],1)./sqrt(size(realPlot_pos_m{ne,ng},1)),'lineprops',{'Color',mmap(np,:)})  
%             title('positive gain')
%         end
%         end
%         xlabel('Time from outcome (s)')
%         ylabel('nanmean z-scored fluorescence')
%         
%         thisPlot = cat(1,realPlot_neg_f{ne,ng}(:,:,[1:end]),realPlot_pos_f{ne,ng}(:,:,[1:end]));
%         thisPlot_m = cat(1,realPlot_neg_m{ne,ng}(:,:,[1:end]),realPlot_pos_m{ne,ng}(:,:,[1:end]));
%         
%         figure(),
%         for np = 1:size(thisPlot,3)
%              subplot(2,2,1); hold on
%             shadedErrorBar(xaxis, nanmean(thisPlot(:,:,np),1),nanstd(thisPlot(:,:,np),[],1)./sqrt(size(thisPlot,1)),'lineprops',{'Color',fmap(np,:)})
%              subplot(2,2,2); hold on
%             shadedErrorBar(xaxis, nanmean(thisPlot_m(:,:,np),1),nanstd(thisPlot_m(:,:,np),[],1)./sqrt(size(thisPlot_m,1)),'lineprops',{'Color',mmap(np,:)})  
%         end
%         xlabel('Time from outcome (s)')
%         ylabel('nanmean z-scored fluorescence')
%         
%         
        
        % Plot combined
%         xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_f{ne,ng},2));
%         thisPlot = cat(1,realPlot_neg_f{ne,ng}(:,:,[1:end]),realPlot_pos_f{ne,ng}(:,:,[1:end]));
%         plotLen = size(realPlot_pos_f{ne,ng},2);
%         nplots = size(thisPlot,3);
%         
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2).*nplots,1);
%         [maxVal,maxTime] = max(thisPlot,[],2);
%         [~,~,whichBin] = histcounts(maxTime,1:plotLen:plotLen*nplots+1);
%         [~,sortIdx] = sort(whichBin);
%         maxVal = maxVal(sortIdx);
%         whichBin = whichBin(sortIdx);
%         sortIdx2 = zeros(size(whichBin));
%         for np = 1:nplots
%             [~,sortIdx2(whichBin==np)] = sort(maxVal(whichBin==np));
%         end
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440 80 531 718]);
%         
%         thisPlot2 = cat(1,predPlot_neg_f{ne,ng}(:,:,[1:end]),predPlot_pos_f{ne,ng}(:,:,[1:end]));
%         thisPlot2 = reshape(thisPlot2,size(thisPlot2,1),size(thisPlot2,2).*nplots,1);
% %         thisPlot = thisPlot(sortIdx,:);
% %         thisPlot2=thisPlot2(sortIdx,:);
% %         sortIdx=sortIdx2;
%         for np = 1:nplots
%             subplot(2,nplots,np); box off
%             imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,(np-1)*plotLen+1:(np)*plotLen),[-1 1]), colormap(flipud(red2blue))
%             ylabel('Neuron (sorted by time of peak activity')
%             title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%             smallcolorbar;
%             xlabel('Time from event')
%              subplot(2,nplots,nplots+np)
%             imagesc(xaxis,1:size(thisPlot2,1),thisPlot2(sortIdx,(np-1)*plotLen+1:(np)*plotLen),[-1 1]), colormap(flipud(red2blue))
%             xlabel('Time from event')
%             smallcolorbar;
%         end
        
      %   print(fullfile(savehere, sprintf('%s_%s_%s_female',savename,events{ne},gains{ne}{ng})), '-dpdf', '-bestfit')

        % Plot male
%         xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_f{ne,ng},2));
%         thisPlot = realPlot_pos_m{ne,ng}(:,:,[1,end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         [~,maxTime] = max(thisPlot,[],2);
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440   378   970   420]);
%         subplot(2,2,1)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_pos_m{ne,ng},2)),[-2.5 2.5]), colormap(flipud(red2blue))
%         ylabel('Neuron (sorted by time of peak activity')
%         title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%         subplot(2,2,2)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_pos_m{ne,ng},2):end),[-2.5 2.5]), colormap(flipud(red2blue))
%         smallcolorbar
%         xlabel('Time from event')
%         thisPlot = predPlot_pos_m{ne,ng}(:,:,[1:end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         subplot(2,2,3)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_pos_m{ne,ng},2)),[-2.5 2.5]), colormap(flipud(red2blue))
%         subplot(2,2,4)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_pos_m{ne,ng},2):end),[-2.5 2.5]), colormap(flipud(red2blue))
%         xlabel('Time from event')
%         smallcolorbar
        
        
        % Plot negative encoding
%          xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_m{ne,ng},2));
%         thisPlot = realPlot_neg_m{ne,ng}(:,:,[1,end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         [~,maxTime] = max(thisPlot,[],2);
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440   378   970   420]);
%         subplot(2,2,1)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_neg_m{ne,ng},2)),[-2.5 2.5]), colormap(flipud(red2blue))
%         ylabel('Neuron (sorted by time of peak activity')
%         title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%         subplot(2,2,2)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_neg_m{ne,ng},2):end),[-2.5 2.5]), colormap(flipud(red2blue))
%         smallcolorbar
%         xlabel('Time from event')
%         thisPlot = predPlot_neg_m{ne,ng}(:,:,[1:end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         subplot(2,2,3)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_neg_m{ne,ng},2)),[-2.5 2.5]), colormap(flipud(red2blue))
%         subplot(2,2,4)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_neg_m{ne,ng},2):end),[-2.5 2.5]), colormap(flipud(red2blue))
%         xlabel('Time from event')
%         smallcolorbar
        
        % Plot combined
%        xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_m{ne,ng},2));
%         thisPlot = cat(1,realPlot_neg_m{ne,ng}(:,:,[1:end]),realPlot_pos_m{ne,ng}(:,:,[1:end]));
%         plotLen = size(realPlot_pos_m{ne,ng},2);
%         nplots = size(thisPlot,3);
%         
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2).*nplots,1);
%         [~,maxTime] = max(thisPlot,[],2);
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440 80 531 718]);
%         
%         thisPlot2 = cat(1,predPlot_neg_m{ne,ng}(:,:,[1:end]),predPlot_pos_m{ne,ng}(:,:,[1:end]));
%         thisPlot2 = reshape(thisPlot2,size(thisPlot2,1),size(thisPlot2,2).*nplots,1);
%            
%         for np = 1:nplots
%             subplot(2,nplots,np)
%             imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,(np-1)*plotLen+1:(np)*plotLen),[-1 1]), colormap(flipud(red2blue))
%             ylabel('Neuron (sorted by time of peak activity')
%             title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%             smallcolorbar;
%             xlabel('Time from event')
%              subplot(2,nplots,nplots+np)
%             imagesc(xaxis,1:size(thisPlot2,1),thisPlot2(sortIdx,(np-1)*plotLen+1:(np)*plotLen),[-1 1]), colormap(flipud(red2blue))
%             xlabel('Time from event')
%             smallcolorbar;
%         end
%            print(fullfile(savehere, sprintf('%s_%s_%s_male',savename,events{ne},gains{ne}{ng})), '-dpdf', '-bestfit')

     
keyboard
    
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
%X = nanzscore(X);
predY = X*cat(1,fit.model_parameters.kernel_profiles{:});

for np = 1:numel(fit.eventNames)
    %extract estimated data by event and gain predictor
    clear thistrials
    %thisEvents = find(fit.kernelPreds{np}==1)+fit.eventDur(np,1)*fit.frameRate;
    thisEvents = fit.event_times{np}';
    for ne = 1:numel(thisEvents)
        if isnan(thisEvents(ne))
            thistrials(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
            thistrials2(ne,:) = nan(size(-fit.eventDur(np,1)*fit.frameRate:fit.eventDur(np,2)*fit.frameRate));
        else
            try
                thistrials(ne,:) = predY(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
                thistrials2(ne,:) = y(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate);
            catch
                thistrials(ne,:) = predY(1,predY(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:end),zeros(thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate-numel(predY),1));
                thistrials2(ne,:) = y(1,y(thisEvents(ne)-fit.eventDur(np,1)*fit.frameRate:end),zeros(thisEvents(ne)+fit.eventDur(np,2)*fit.frameRate-numel(y),1));
                
            end
        end
    end
%     if np==5
%         keyboard
%     end
    predY_event{np} = nanmean(thistrials);
    semPredY_event{np} = nanstd(thistrials)./sqrt(sum(~isnan(thisEvents)));
    for ng = 1:numel(fit.gainNames{np})
        [~,~,bins] = histcounts(fit.gain_preds{np,ng},fit.gainBins{np,ng});
       % bins = bins(~isnan(fit.event_times{np}));
        for nb = 1:numel(fit.gainBins{np,ng})
            predFluor_bins{np,nb,ng} = nanmean(thistrials(bins==nb,:),1);
            predFluor_value{np,nb,ng} = nanmean(thistrials(fit.gain_preds{np,ng}==fit.gainBins{ng}(nb),:),1);
            semPredFluor_bins{np,nb,ng} = nanstd(thistrials(bins==nb,:),[],1)./sqrt(sum(bins==nb));
            semPredFluor_value{np,nb,ng} = nanstd(thistrials(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb),:),[],1)./sqrt(sum(fit.gain_preds{np,ng}==fit.gainBins{np,ng}(nb)));

        end
    end
end


%         plotLen = size(realPlot_pos_f{ne,ng},2);
%         nplots = size(thisPlot,3);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*nplots,1);
%         [~,maxTime] = max(thisPlot,[],2);
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440   378   970   420]);
%         
%         for np=1:nplots    
%         subplot(2,nplots,np)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_pos_f{ne,ng},2)),[-1.5 1.5]), colormap(flipud(red2blue))
%         ylabel('Neuron (sorted by time of peak activity')
%         title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%         subplot(2,2,2)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_pos_f{ne,ng},2):end),[-1.5 1.5]), colormap(flipud(red2blue))
%         smallcolorbar
%         xlabel('Time from event')
%         thisPlot = predPlot_pos_f{ne,ng}(:,:,[1 end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         subplot(2,2,3)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_pos_f{ne,ng},2)),[-1.5 1.5]), colormap(flipud(red2blue))
%         subplot(2,2,4)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_pos_f{ne,ng},2):end),[-1.5 1.5]), colormap(flipud(red2blue))
%         xlabel('Time from event')
%         smallcolorbar
%         end
%         
%         % Plot negative encoding
%          xaxis = linspace(-eventDur(ne,1), eventDur(ne,2),size(realPlot_neg_f{ne,ng},2));
%         thisPlot = realPlot_neg_f{ne,ng}(:,:,[1,end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         [~,maxTime] = max(thisPlot,[],2);
%         [~,sortIdx] = sort(maxTime,'ascend');
%         figure('Position',[440   378   970   420]);
%         subplot(2,2,1)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_neg_f{ne,ng},2)),[-1.5 1.5]), colormap(flipud(red2blue))
%         ylabel('Neuron (sorted by time of peak activity')
%         title(sprintf('event: %s gain: %s',events{ne},gains{ne}{ng}))
%         subplot(2,2,2)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_neg_f{ne,ng},2):end),[-1.5 1.5]), colormap(flipud(red2blue))
%         smallcolorbar
%         xlabel('Time from event')
%         thisPlot = predPlot_neg_f{ne,ng}(:,:,[1 end]);
%         thisPlot = reshape(thisPlot,size(thisPlot,1),size(thisPlot,2)*2,1);
%         subplot(2,2,3)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1:size(realPlot_neg_f{ne,ng},2)),[-1.5 1.5]), colormap(flipud(red2blue))
%         subplot(2,2,4)
%         imagesc(xaxis,1:size(thisPlot,1),thisPlot(sortIdx,1+size(realPlot_neg_f{ne,ng},2):end),[-1.5 1.5]), colormap(flipud(red2blue))
%         xlabel('Time from event')
%         smallcolorbar
end