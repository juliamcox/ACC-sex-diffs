function bilinearRegression_summary(regVer,savehere,rawFlag)

% Parameters
fbasename = fullfile(whereAreWe('bilinear'),savehere);
% Regression parameters 
[events, kernelType, gains, eventDur,nbasis] = getEvents_bilinear(regVer);
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
sigLevel  = .01;
try load(fullfile(fbasename,regVer,'summary.mat'))
catch
    % load male and female animal IDs
    load(fullfile(fbasename,'recs'));
    load(fullfile(whereAreWe('imaging'),'activeNeurons_f.mat'));
    activeNeurons_f = activeNeurons;
    load(fullfile(whereAreWe('imaging'),'activeNeurons_m.mat'));
    activeNeurons_m = activeNeurons;
    
    
    % Filenames
    cd(fullfile(fbasename,regVer));
    fnames = dir;
    fnames = {fnames(:).name};
    fnames = fnames(contains(fnames,'roi'));
    fnames = fnames(contains(fnames,sprintf('rawFlag%d',rawFlag)));
    % sort fnames by recs and roi number so that same order as activeNeurons
    thisFnames = [];
    for na = 1:numel(recs_f)
        temp = fnames(contains(fnames,recs_f{na}));
        idx1 = strfind(temp,'roi');
        idx2 = strfind(temp,'_');
        roi = cell2mat(cellfun(@(x, y, z) str2num(x(y+3:z(1)-1)), temp,idx1,idx2,'UniformOutput',false));
        [~,idx3] = sort(roi);
        thisFnames = cat(2,thisFnames,temp(idx3)); % find
    end
    thisFnames_m = [];
    for na = 1:numel(recs_m)
        temp = fnames(contains(fnames,recs_m{na}));
        idx1 = strfind(temp,'roi');
        idx2 = strfind(temp,'_');
        roi = cell2mat(cellfun(@(x, y, z) str2num(x(y+3:z(1)-1)), temp,idx1,idx2,'UniformOutput',false));
        [~,idx3] = sort(roi);
        thisFnames_m = cat(2,thisFnames_m,temp(idx3)); % find
    end
    
    % Extract data
    try load(fullfile(fbasename,regVer,'summary.mat'),'summary_f','summary_m')
    catch
        summary_f = extractData(recs_f,thisFnames,savehere,regVer);
        summary_m = extractData(recs_m,thisFnames_m,savehere,regVer);
        save(fullfile(fbasename,regVer,'summary.mat'),'summary_f','summary_m');
    end
    % Plot proportion siginifcant for each event and gain
    model_pvals = cell2mat(summary_f.model_pval');
    model_r2 = cell2mat(summary_f.model_r2');
    model_r2_orig = cell2mat(summary_f.model_r2_orig');
    model_pvals_m = cell2mat(summary_m.model_pval');
    model_r2_m = cell2mat(summary_m.model_r2');
    model_r2_orig_m = cell2mat(summary_m.model_r2_orig');
    clear pvals pvals_m
    idx = 1;
    for ne = 1:numel(events)
        for ng = 1:numel(gains{ne})+1
            pvals(:,idx) = cell2mat(cellfun(@(x) x(:,ng,ne), summary_f.gain_pval,'UniformOutput',false));
            pvals_m(:,idx) = cell2mat(cellfun(@(x) x(:,ng,ne), summary_m.gain_pval,'UniformOutput',false));
            
            % extract coefficients
            coeffs{ng,ne} = cell2mat(cellfun(@(x) x(:,ng,ne), summary_f.gain_coeff,'UniformOutput',false));
            coeffs_m{ng,ne} = cell2mat(cellfun(@(x) x(:,ng,ne), summary_m.gain_coeff,'UniformOutput',false));
            idx = idx+1;
        end
    end
    try
        pvals = pvals(activeNeurons_f==1,:);
        pvals_m = pvals_m(activeNeurons_m==1,:);
    catch
        keyboard
    end
    idx = 1;
    for ne = 1:numel(events)
        for ng = 1:numel(gains{ne})+1
            numSig(ng,ne) = sum(pvals(:,idx)<sigLevel);
            numNotSig(ng,ne) =sum(pvals(:,idx)>=sigLevel);
            numSig_m(ng,ne) = sum(pvals_m(:,idx)<sigLevel);
            numNotSig_m(ng,ne) =sum(pvals_m(:,idx)>=sigLevel);
            idx =idx+1;
        end
    end
    
    
    
    % Test for significant differences in proportion significantly encoding each event
    for ne = 1:numel(events)
        for ng = 1:numel(gains{ne})+1
            [chi2_pval(ne,ng),chi2stat(ne,ng)] = chi2PropTest(numSig(ng,ne),numSig_m(ng,ne),size(pvals,1),size(pvals_m,1),0);
        end
    end
    
    thisFnames = thisFnames(boolean(activeNeurons_f));
    thisFnames_m = thisFnames_m(boolean(activeNeurons_m));
    save(fullfile(fbasename,regVer,'summary.mat'),'fnames','thisFnames','thisFnames_m','chi2_pval','chi2stat','numSig','numSig_m','numNotSig','numNotSig_m','pvals','pvals_m','model_r2','model_r2_m','coeffs','coeffs_m','-append');
end
totGains = size(coeffs,1);


f=figure('Position',get(0,'ScreenSize'));
for ne = 1:numel(events)
    for ng = 1:numel(gains{ne})+1
        if ng<=numel(gains{ne})
            thisplot = ng;
        else
            thisplot = size(coeffs,1);
        end
        subplot(totGains,numel(events),ne+(thisplot-1)*numel(events)); hold on
        b=bar(1,numSig(ng,ne)./(numSig(ng,ne)+numNotSig(ng,ne)));
        hold on, b(2) = bar(2,numSig_m(ng,ne)./(numSig_m(ng,ne)+numNotSig_m(ng,ne)));
        b(1).FaceColor = plotParams.femaleC;
        b(1).EdgeColor = 'none';
        b(2).FaceColor = plotParams.maleC;
        b(2).EdgeColor = 'none';
        if chi2_pval(ne,ng)<sigLevel
            title(sprintf('%s^*',events{ne}))
        else
            title(events{ne})
        end
        g = gca;
        try
            ylabel(gains{ne}{ng})
            g.YLim = [0 .16];
        catch
            ylabel('offset')
        end
        
        g.XAxis.Color = 'none';
        %g.YAxis.Color = 'none';
        g.YAxis.Label.Color = [0 0 0];
        g.YAxis.Label.Visible = 'on';
        if ng<totGains
            %  g.YLim = [0 .16];
        end
    end
end




function summary = extractData(recs,fnames,savehere,regVer)


% base filename
fbasename = fullfile(whereAreWe('imaging'),'bilinearRegression',savehere,regVer);
% load regression parameters
[events, kernelType, gains, eventDur,nbasis] = getEvents_bilinear(regVer);
% initialize output structure

summary.gain_pval     = cell(size(recs));
summary.model_pval    = cell(size(recs));
summary.model_r2      = cell(size(recs));
summary.gain_coeff    = cell(size(recs));
summary.model_r2_orig = cell(size(recs));


for na = 1:numel(recs)
    thisFnames = fnames(contains(fnames,recs{na})); % find rois for each recording
    summary.model_pval{na} = zeros(size(thisFnames));
    summary.model_r2{na} = zeros(size(thisFnames));
    
    for nf = 1:numel(thisFnames)
        load(fullfile(fbasename,thisFnames{nf}),'fit_crossVal','pvalues','pvalues_coeff','pvalues_crossVal','pvalues_coeff_crossVal')
        summary.model_r2{na}(nf) = mean(cell2mat(cellfun(@(x) x.R2(end), fit_crossVal.T, 'UniformOutput', false)));
        summary.model_pval{na}(nf) = mean(cell2mat(cellfun(@(x) x.pValue(end), fit_crossVal.T, 'UniformOutput', false)));
        summary.model_r2_orig{na}(nf) = mean(fit_crossVal.orig.T.R2);
        
        for ne = 1:numel(events)
            summary.gain_pval{na}(nf,:,ne) = pvalues(ne,:);
            thisCoeff = cell2mat(cellfun(@(x) x.trial_coefficients{ne},fit_crossVal.model_parameters,'UniformOutput',false));
            thisCoeff = mean(thisCoeff,2);
            summary.gain_coeff{na}(nf,1:numel(thisCoeff),ne) = thisCoeff;
            summary.gain_pval_crossVal{na}(nf,:,ne) = pvalues_crossVal(ne,:);
            summary.pval_coeff_crossVal{na}(nf,:,ne) = pvalues_coeff_crossVal(ne,:);
            summary.pval_coeff{na}(nf,:,ne) = pvalues_coeff(ne,:);
        end
        
    end
end