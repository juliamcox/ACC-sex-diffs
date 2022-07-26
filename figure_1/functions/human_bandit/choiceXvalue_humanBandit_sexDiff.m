function choiceXvalue_humanBandit_sexDiff(cohort,nbins,perfThresh,sideThresh)

% cohort: which subjects to plot
% qLoc: location of q-learning stan fit 



basefilename = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit'); 

%% Load each subject's data extract value bins 
sList_f = dir(fullfile(basefilename, [cohort '_f'], 'qLearn*.mat')); 
sList_f = {sList_f(:).name};
sList_m = dir(fullfile(basefilename, [cohort '_m'], 'qLearn*.mat')); 
sList_m = {sList_m(:).name};


[choice_f,bins]=extractChoiceXValue(sList_f,[cohort '_f'],nbins,perfThresh,sideThresh);
save(fullfile(basefilename,cohort, sprintf('choiceXvalue_f_bins%d_perfThresh%s_sideThresh%s.mat',nbins,num2str(perfThresh),num2str(sideThresh))),'choice_f','bins');
[choice_m,bins]=extractChoiceXValue(sList_m,[cohort '_m'],nbins,perfThresh,sideThresh);
save(fullfile(basefilename,cohort, sprintf('choiceXvalue_m_bins%d_perfThresh%s_sideThresh%s.mat',nbins,num2str(perfThresh),num2str(sideThresh))),'choice_m','bins');




%% 
function [choice, bins] = extractChoiceXValue(fList,cohort,nbins,perfThresh,sideThresh)


valTypes = {'QDiff';'QChosenDiff';'QTot';'QChosen';'QUnchosen';'QRight';'QLeft'}; 

bins.QDiff = linspace(-1,1,nbins);
bins.QChosenDiff = bins.QDiff;
bins.QTot = linspace(0,2,nbins);
bins.QChosen = linspace(0,1,nbins);
bins.QUnchosen = bins.QChosen;
bins.QRight = bins.QChosen;
bins.QLeft = bins.QChosen;


basefilename = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit',cohort); 

x = 1;
for ns = 1:numel(fList)
    load(fullfile(basefilename,fList{ns})); % load q-values
    load(fullfile(basefilename,fList{ns}(strfind(fList{ns},'_')+1:end))); % load behavior
    
    thisStay = [NaN; data.choice(1:end-1)==data.choice(2:end)];
    prevRew = [NaN; data.reward(1:end-1)];
    if mean(thisStay(prevRew==1))-mean(thisStay(prevRew==0)) >= perfThresh && abs(mean(qLearn.side))<sideThresh && numel(data.choice)>=200
        try
            choice.age(x)   = str2num(data.raw{end}.responses.age);
        catch
            choice.age(x) = NaN;
        end
        choice.aID(x)   = ns;
        for nv = 1:numel(valTypes)
            thisvalue = eval(sprintf('qLearn.%s',valTypes{nv}));
            % bins for quantiles
            tempBins = prctile(thisvalue,[0 100/nbins:100/nbins:100]);
            [~,~,quantBins] = histcounts(thisvalue,nbins);
            
            thisbins  = eval(sprintf('bins.%s',valTypes{nv}));
            [~,~,thisvalueBins] = histcounts(thisvalue,thisbins);
            for nb = 1:numel(thisbins)
                eval(sprintf('choice.probRight_%s(x,nb) = nanmean(data.choice(thisvalueBins==nb)==0);',valTypes{nv}));
                eval(sprintf('choice.probRight_%s_quant(x,nb) = nanmean(data.choice(quantBins==nb)==0);',valTypes{nv}));
                
            end
            
            eval(sprintf('choice.%s_trials{x} = thisvalue;',valTypes{nv}));
            
        end
        choice.choice_trials{x} = data.choice;
        x = x+1;
    end
end