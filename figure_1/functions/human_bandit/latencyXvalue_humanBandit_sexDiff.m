function latencyXvalue_humanBandit_sexDiff(cohort,nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh,sList_f,sList_m)

% cohort: which subjects to plot
% qLoc: location of q-learning stan fit 


basefilename = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit'); 

%% Load each subject's data extract value bins 
if nargin < 9
    sList_f = dir(fullfile(basefilename, [cohort '_f'], ['qlearn*.mat']));
    sList_f = {sList_f(:).name};
    sList_m = dir(fullfile(basefilename, [cohort '_m'], ['qlearn*.mat']));
    sList_m = {sList_m(:).name};
end


[latency_f,bins,alpha,beta,stay,side]=extractLatencyXValue(sList_f,[cohort '_f'],nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
save(fullfile(basefilename,cohort, sprintf('latencyXvalue_f_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_f','bins','alpha','beta','side','stay');
[latency_m,bins,alpha_m,beta_m,stay_m,side_m]=extractLatencyXValue(sList_m,[cohort '_m'],nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh);
save(fullfile(basefilename,cohort, sprintf('latencyXvalue_m_bins%d_lowCut%d_highCut%d_zscore%d_perfThresh%d_sideThresh%d.mat',nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh*10,sideThresh*10)),'latency_m','bins','alpha_m','beta_m','side_m','stay_m');


end
%% 
function [latency, bins, alpha,beta,stay,side] = extractLatencyXValue(fList,cohort,nbins,lowCutoff,highCutoff,zscoreFlag,perfThresh,sideThresh)


valTypes = {'QDiff';'QChosenDiff';'QTot';'QChosen';'QUnchosen'}; 
latencyTypes = {'choice';'trialStart'}; 

bins.QDiff = linspace(-1,1,nbins+1);
bins.QChosenDiff = bins.QDiff;
bins.QTot = linspace(0,2,nbins+1);
bins.QChosen = linspace(0,1,nbins+1);
bins.QUnchosen = bins.QChosen;

% create pathname for q-values 
basefilename = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit',cohort); 
latency.aList = [];

acounter = 1; 
for ns = 1:numel(fList)
    load(fullfile(basefilename,fList{ns})); % load q-values
    load(fullfile(basefilename,fList{ns}(strfind(fList{ns},'_')+1:end))); % load behavior
    
    thisStay = [NaN; data.choice(1:end-1)==data.choice(2:end)];
    prevRew = [NaN; data.reward(1:end-1)];
    if mean(thisStay(prevRew==1))-mean(thisStay(prevRew==0)) >= perfThresh && abs(mean(qLearn.side))<sideThresh && numel(data.choice)>200
        try
            latency.age(acounter)   = str2num(data.raw{end}.responses.age);
        catch
            latency.age(acounter) = NaN;
            keyboard
        end
        latency.aList = cat(1,latency.aList,fList(ns));
        alpha(acounter) = mean(qLearn.alpha);
        beta(acounter)  = mean(qLearn.beta);
        side(acounter)  = mean(qLearn.side);
        stay(acounter)  = mean(qLearn.stay);
        
        for nl = 1:numel(latencyTypes)
            
            thislatency = eval(sprintf('data.%sLatency',latencyTypes{nl}));
            thislatency(thislatency<lowCutoff | thislatency > highCutoff) = NaN;
            if zscoreFlag == 1
                thislatency = nanzscore(thislatency);
            elseif zscoreFlag == 2
                dim = find(max(size(thislatency)));
                mu = nanmean(thislatency,dim);
                thislatency = (thislatency - mu);
            end
            
            %thislatency(1:20) = NaN;
            
            if isempty(thislatency)
                for nv = 1:numel(valTypes)
                    thisvalue = eval(sprintf('qLearn.%s',valTypes{nv}));
                    
                    % bins for quantiles
                    tempBins = prctile(thisvalue,[0 100/binNum:100/binNum:100]);
                    [~,~,quantBins] = histcounts(thisvalue,tempBins);
                    % bins for value
                    thisbins  = eval(sprintf('bins.%s',valTypes{nv}));
                    [~,~,thisvalueBins] = histcounts(thisvalue,thisbins);
                    for nb = 1:numel(thisbins)
                        eval(sprintf('latency.%s_%s(acounter,nb) = NaN;',latencyTypes{nl},valTypes{nv}));
                        eval(sprintf('latency.%s_%s_quant(acounter,nb) = NaN;',latencyTypes{nl},valTypes{nv}));
                        %eval(sprintf('latency.%s_%s_trials = {acounter,nb} = NaN;' ,latencyTypes{nl},valTypes{nv}));
                    end
                    latency.valueBins = thisvalueBins;
                    latency.quantBins = quantBins;
                    latency.trials = eval(latencyTypes{nl});
                    
                end
            else
                for nv = 1:numel(valTypes)
                    thisvalue = eval(sprintf('qLearn.%s',valTypes{nv}));
                    % bins for quantiles
                    tempBins = prctile(thisvalue,[0 100/nbins:100/nbins:100]);
                    [~,~,quantBins] = histcounts(thisvalue,tempBins);
                    % bins for value
                    thisbins  = eval(sprintf('bins.%s',valTypes{nv}));
                    [~,~,thisvalueBins] = histcounts(thisvalue,thisbins);
                    for nb = 1:numel(thisbins)-1
                        eval(sprintf('latency.%s_%s(acounter,nb) = nanmean(thislatency(thisvalueBins==nb));',latencyTypes{nl},valTypes{nv}));
                        eval(sprintf('latency.%s_%s_quant(acounter,nb) = nanmean(thislatency(quantBins==nb));',latencyTypes{nl},valTypes{nv}));
                        eval(sprintf('latency.%s_%s_count(acounter,nb) = sum(thisvalueBins==nb);',latencyTypes{nl},valTypes{nv}));
                        eval(sprintf('latency.%s_%s_trials{acounter,nb} = thislatency(thisvalueBins==nb);',latencyTypes{nl},valTypes{nv}));
                    end
                    eval(sprintf('latency.%s_value{acounter} = thisvalue;',valTypes{nv}));
                    eval(sprintf('latency.%s_valueQuant{acounter} = quantBins;',valTypes{nv}));
                end
            end
            eval(sprintf('latency.%s_trials{acounter} = thislatency;',latencyTypes{nl}));
        end
        acounter = acounter+1;
    end
    
    
end
end