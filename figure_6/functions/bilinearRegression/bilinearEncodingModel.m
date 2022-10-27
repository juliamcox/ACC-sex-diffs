function fit = bilinearEncodingModel(recIdx,savehere,regVer,plotFlag,rawFlag,frameRate,qFile)
% recIdx = index of neuron 
% savehere = location to save fit and load data files
% regVer = which version of the regression, see: getEvents_bilinear
% plotFlag = plot estimated and real data
% rawFlag = use raw or denoised output from CNMFe
% frameRate = frame rate of the recording
% qFile = q-model version 
tic

%% Parameters

fbasename = fullfile(whereAreWe('bilinear'),savehere,sprintf('rawFlag%d_frameRate%d_%s',rawFlag,frameRate,qFile(1:end-4)));
load(fullfile(fbasename,'fnames.mat'));
rec = fnames{recIdx};
fprintf('%s \n',rec)
underIdx  = strfind(rec,'_');
nCVFolds  = 5;
numShuff  = 1000;


% load parameters for regression
[events, kernelType, gains, eventDur,nbasis,numConsec] = getEvents_bilinear(regVer); 

% create output structure with model parameters
fit.kernelType = kernelType;
fit.frameRate  = frameRate;
fit.rawFlag    = rawFlag;
fit.eventNames = events;
fit.gainNames  = gains;
fit.eventDur   = eventDur;
fit.numConsec  = numConsec;
if contains(kernelType,'timeLag')
    for ne = 1:numel(events)
        nbasis(ne) = numel(-eventDur(ne,1):eventDur(ne,2));
    end
    fit.nbasis = nbasis; 
else
    fit.nbasis = nbasis;
end

%% Load data
y=load(fullfile(fbasename,rec),'y'); % load fluorescence
y = y.y;
y = nanzscore(y);
recLoc = fullfile(whereAreWe('imaging'),rec(underIdx(1)+1:underIdx(2)-1),rec(underIdx(2)+1:underIdx(5)-1));
load(fullfile(recLoc,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']));
valueFname = fullfile(recLoc,sprintf('valueData_fs%d_raw%d_nfolds%d_%s',frameRate,rawFlag,nCVFolds,qFile));


%% Generate predictors 
fit = extractPredictors(fit,tdtEvents,valueFname);

%% Extract event-triggered fluorescence 
fit = extractEvents_bilinear(fit,y); 

%% Fit model with cross validation to calculate goodness of fit 
fit_crossVal = bilinearEncodingModel_fit_crossVal(fit,y);

%% Fit model to all data and shifted data to determine significance of coefficients 

%%% Fit full model 
fit_real = bilinearEncodingModel_fit(fit,y);

%% Fit on shuffled data
% 
p = gcp('nocreate'); % If no pool, do not create new one.
delete(p);
if isempty(p)
    try
        pc = parcluster('local');
        % delete(pc.Jobs);
        if contains(pc.Host,'scotty')
            pc.NumWorkers = 44;
            parpool(pc,44);
        else
            pc.NumWorkers = 24;
            parpool(pc,24);
        end
    catch
        parpool(12);
    end
else
    delete(gcp('nocreate'))
    try
        pc = parcluster('local');
        pc.NumWorkers = 24;
        
        parpool(pc,24);
     catch
         parpool(12);
     end
end
fit_shuff = cell(numShuff,1);
parfor ns = 1:numShuff
    pause(rand(1))
    thisFit = fit;
    thisy = circshift(y,randi(numel(y)-thisFit.frameRate,1));
    fit_shuff{ns} = bilinearEncodingModel_fit(thisFit,thisy);
fit_shuff{ns}.y_pred = [];
fit_shuff{ns}.X = [];
end
temp.tstat = cellfun(@(x) x.tstat,fit_shuff,'UniformOutput',false);
temp.model_parameters = cellfun(@(x) x.model_parameters,fit_shuff,'UniformOutput',false);
fit_shuff = temp;
%% Calculate p-values
% 
for ne = 1:numel(fit.eventNames)
    for ng = 1:numel(fit.gainNames{ne})+1
        thisT_shuff = cellfun(@(x) x{ne}(ng),fit_shuff.tstat);   
        thisT_real = fit_real.tstat{ne}(ng);
        if thisT_real < 0 
            pvalues(ne,ng) = 1-mean(thisT_real<thisT_shuff);
        else
            pvalues(ne,ng) = 1-mean(thisT_real>thisT_shuff);
        end
        
        thisCoeff_real = fit_real.model_parameters.trial_coefficients{ne}(ng);
        thisCoeff_shuff = cellfun(@(x) x.trial_coefficients{ne}(ng),fit_shuff.model_parameters);
        if thisT_real < 0
            pvalues_coeff(ne,ng) = 1-mean(thisCoeff_real<thisCoeff_shuff);
        else
            pvalues_coeff(ne,ng) = 1-mean(thisCoeff_real>thisCoeff_shuff);
        end
        
        
        thisT_shuff = cellfun(@(x) x{ne}(ng),fit_shuff.tstat);
        thisT_real = mean(cellfun(@(x) x{ne}(ng), fit_crossVal.tstat));
        if thisT_real < 0
            pvalues_crossVal(ne,ng) = 1-mean(thisT_real<thisT_shuff);
        else
            pvalues_crossVal(ne,ng) = 1-mean(thisT_real>thisT_shuff);
        end
        
        thisCoeff_real = mean(cellfun(@(x) x.trial_coefficients{ne}(ng),fit_crossVal.model_parameters));
        thisCoeff_shuff = cellfun(@(x) x.trial_coefficients{ne}(ng),fit_shuff.model_parameters);
        if thisT_real < 0
            pvalues_coeff_crossVal(ne,ng) = 1-mean(thisCoeff_real<thisCoeff_shuff);
        else
            pvalues_coeff_crossVal(ne,ng) = 1-mean(thisCoeff_real>thisCoeff_shuff);
        end
        
    end
end


%% Find mean cross-validated r2
for nr = 1:numel(fit_crossVal.T)
   r2(nr) = fit_crossVal.T{nr}.R2(end);
end
r2=mean(r2);

%% Plot
if plotFlag
   plotBilinearRegression(fit_real,r2); 
end

%% save data
if ~isdir(fullfile(whereAreWe('imaging'),'bilinearRegression',savehere,regVer))
    mkdir(fullfile(whereAreWe('imaging'),'bilinearRegression',savehere,regVer))
end
save(fullfile(whereAreWe('imaging'),'bilinearRegression',savehere,regVer,rec),'fit_crossVal','fit_real','fit_shuff','pvalues','pvalues_coeff','pvalues_coeff_crossVal','pvalues_crossVal','-v7.3');
toc


delete(pc);
 %% trashcan
% %%% Fit reduced model (events)
% for ne = 1:numel(fit.eventNames)
%     % Generate data structure for reduced model 
%     fit_reduced.kernelType = kernelType;
%     fit_reduced.frameRate  = frameRate;
%     fit_reduced.rawFlag    = rawFlag;
%     fit_reduced.eventNames = events(~contains(events,fit.eventNames{ne}));
%     fit_reduced.gainNames  = gains;
%     fit_reduced.eventDur   = eventDur(~contains(events,fit.eventNames{ne}),:);
%     if contains(kernelType,'timeLag')
%         for ne = 1:numel(events)
%             nbasis(ne) = numel(-eventDur(ne,1):eventDur(ne,2));
%         end
%         fit_reduced.nbasis = nbasis;
%     else
%         fit_reduced.nbasis = nbasis(~contains(events,fit.eventNames{ne}));
%     end
%     fit_reduced = extractPredictors(fit_reduced,tdtEvents,valueFname);
%     fit_reduced = bilinearEncodingModel_fit(fit_reduced,y); 
%     % Calculate f-stat
%     fstat_event(ne) = getFstat_bilinear(fit_full,fit_reduced,y);
% end
% 
% %%% Fit reduced model (gains)
% for ne = 1:numel(fit.eventNames)
%     for ng = 1:numel(fit.gainNames)
%         % Generate data structure for reduced model
%         fit_reduced.kernelType = kernelType;
%         fit_reduced.frameRate  = frameRate;
%         fit_reduced.rawFlag    = rawFlag;
%         fit_reduced.eventNames = events;
%         fit_reduced.gainNames  = gains{ne}(setdiff(1:numel(gains{ne}),ng));
%         fit_reduced.eventDur   = eventDur(~contains(events,fit.eventNames{ne}),:);
%         if contains(kernelType,'timeLag')
%             for ne = 1:numel(events)
%                 nbasis(ne) = numel(-eventDur(ne,1):eventDur(ne,2));
%             end
%             fit_reduced.nbasis = nbasis;
%         else
%             fit_reduced.nbasis = nbasis(~contains(events,fit.eventNames{ne}));
%         end
%         fit_reduced = extractPredictors(fit_reduced,tdtEvents,valueFname);
%         fit_reduced = bilinearEncodingModel_fit(fit_reduced,y);
%         % Calculate f-stat
%         fstat_gain(ne,ng) = getFstat_bilinear(fit_full,fit_reduced,y);
%     end
% end
% 
% %%% Fit reduced model (events)
%     for ne = 1:numel(fit.eventNames)
%         % Generate data structure for reduced model
%         
%         fit_reduced.eventNames = events(~contains(events,fit.eventNames{ne}));
%         fit_reduced.gainNames  = gains;
%         fit_reduced.eventDur   = eventDur(~contains(events,fit.eventNames{ne}),:);
%         
%         fit_reduced = extractPredictors(fit_reduced,tdtEvents,valueFname);
%         fit_reduced = bilinearEncodingModel_fit(fit_reduced,y);
%         % Calculate f-stat
%         fstat_event_shuffle{ns}(ne) = getFstat_bilinear(fit_shuff{ns},fit_reduced,y);
%     end
%     
%     %%% Fit reduced model (gains)
%     for ne = 1:numel(fit.eventNames)
%         for ng = 1:numel(fit.gainNames)
%             % Generate data structure for reduced model
%             fit_reduced.kernelType = kernelType;
%             fit_reduced.frameRate  = frameRate;
%             fit_reduced.rawFlag    = rawFlag;
%             fit_reduced.eventNames = events;
%             fit_reduced.gainNames  = gains{ne}(setdiff(1:numel(gains{ne}),ng));
%             fit_reduced.eventDur   = eventDur(~contains(events,fit.eventNames{ne}),:);
%             
%             fit_reduced = extractPredictors(fit_reduced,tdtEvents,valueFname);
%             fit_reduced = bilinearEncodingModel_fit(fit_reduced,y);
%             % Calculate f-stat
%             fstat_gain_shuffle{ns}(ne,ng) = getFstat_bilinear(fit_shuff{ns},fit_reduced,y);
%         end
%     end