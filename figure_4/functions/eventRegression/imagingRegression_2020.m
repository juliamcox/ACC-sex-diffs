function imagingRegression_2020(rec,regName,zscoreFlag,rawFlag,frameRate,statMethod,num_shuff,reStart,expt)
%% Figure 4 



% Run this function for each recording 

% Event-based, regression a la Parker et al 2020 inspired by ben_regression_matform.m; ben_regression_script; get_f_pvals_lasso.m; get_f_pvals_reg.m 
% 25 degree of freedom spline basis set that spanned -2 to 6 sec for action events and 0 to 8s for stimulus events


%rec: file name/location
%regName: which predictors? see getEvents.m
%zscoreFlag: zscore gcamp data and predictor matrix?
%rawFlag: use raw or denoised gcamp data? if rawFlag == 2 use deconvolved spike times
%frameRate: frame rate for regression 
%statMethod = how to determine significant events? %f (calculate f-statistic between full and reduced model w/o predictors for a given event); crossval (compares goodness of fit (R^2) between full and reduced models across runs to determine significance) 
%num_shuff: number of shuffles for null model 
%restart: start from loaded data or from beginning
%% Parameters
if nargin < 8
    expt = 'ACC_DMS_imaging';
elseif nargin < 7
    reStart = 0;
elseif nargin < 6
    num_shuff = 500;
    reStart = 0; 
end

fext = 'fig4';
nCVFolds = 5; 
fext = ['basisSet_' fext '_' regName];

%%% Set pathname
fbasename = fullfile(whereAreWe('imaging'), expt); 



%% Extract metadata
info = matfile(fullfile(fbasename, rec, 'info.mat'));
frameRateO = info.frameRate;
droppedIdx = info.droppedIdx;
try
    nframesLog = info.nframesLog;
catch
    nframesLog = info.nFramesLog;
end

%filename extension for saving
fnext = sprintf('%s_fs%d_raw%d_zscore%d_%s','linReg', frameRate, rawFlag,zscoreFlag,fext);




%% Load and downsample gcamp signal as needed
%load gcamp data
if rawFlag == 2
    try dff = matfile(fullfile(fbasename,rec,'analyzed_neuron_final.mat'));
        dff = eval(sprintf('dff.dff_s_%d',frameRate));
    catch
        load(fullfile(fbasename,rec,'analyzed_neuron_final.mat'),'neuron');
        dff = neuron.S;
        if size(dff,2) > size(dff,1)
            dff = dff';
        end
        %if droppped frames had been removed for cnmfe, add back in NaN
        if numel(droppedIdx)>0 numel(droppedIdx)+nframesLog > length(dff)
            temp = zeros(numel(droppedIdx)+nframesLog,size(dff,2));
            temp(droppedIdx,:) = NaN;
            temp(setdiff(1:numel(droppedIdx)+nframesLog,droppedIdx),:) = dff;
            dff = temp;
        elseif numel(droppedIdx) == 0 && nframesLog ~= length(dff)
            keyboard %this is probably T55 7/24 where I had to cut two frames to get inscopix to export tiff
        end
        %downsample if necessary
        if frameRate~=frameRateO
            dff = downsampleImaging(frameRateO, frameRate,dff);
        end
        eval(sprintf('dff_s_%d = dff;', frameRate));
        save(fullfile(fbasename,rec,'analyzed_neuron_final.mat'), ['dff_s_' num2str(frameRate)], '-append');
    end
else
    try dff = matfile(fullfile(fbasename,rec,'analyzed_neuron_final.mat'));
        dff = eval(sprintf('dff.dff_%d_raw%d',frameRate,rawFlag));
    catch
        load(fullfile(fbasename,rec,'analyzed_neuron_final.mat'),'neuron');
        if rawFlag %if rawFlag use the raw fluorescence
            dff = neuron.C_raw;
        else %otherwise use denoised fluorescence
            dff = neuron.C;
        end
        
        if size(dff,2) > size(dff,1)
            dff = dff';
        end
        %if droppped frames had been removed for cnmfe, add back in NaN
        if droppedIdx == 0
            droppedIdx = [];
        end
        if numel(droppedIdx)>0 numel(droppedIdx)+nframesLog > length(dff)
            temp = zeros(numel(droppedIdx)+nframesLog,size(dff,2));
            temp(droppedIdx,:) = NaN;
            temp(setdiff(1:numel(droppedIdx)+nframesLog,droppedIdx),:) = dff;
            dff = temp;
        elseif numel(droppedIdx) == 0 && nframesLog ~= length(dff)
            keyboard %this is probably T55 7/24 where I had to cut two frames to get inscopix to export tiff
        end
        %downsample if necessary
        if frameRate~=frameRateO
            fprintf('Downsampling....\n')
            dff = downsampleImaging(frameRateO, frameRate,dff);
        end
        eval(sprintf('dff_%d_raw%d = dff;', frameRate,rawFlag));
        save(fullfile(fbasename,rec,'analyzed_neuron_final.mat'), sprintf('dff_%d_raw%d', frameRate,rawFlag), '-append');
    end
end

%% Generate predictor matrix
nframes = size(dff,1);
try load(fullfile(fbasename, rec, ['M_' fnext '.mat']));
   X = x_basic;
catch
    [X,con_iden] = createPredictorMatrix_2020(rec, fnext,frameRate,frameRateO,regName,nframes,expt);
end


%get indices for events
num_cons = unique(con_iden);
for n = 1:numel(num_cons)
    preds_cells{n} = find(con_iden==n);
end

% normalize gcamp and predictor matrix
if zscoreFlag %lasso automatically standardizes, although coefficients output in original scale
    mu = nanmean(dff,1);
    sigma = nanstd(dff,[],1);
    gcamp = (dff-repmat(mu,size(dff,1),1))./repmat(sigma, size(dff,1),1);
    X = zscore(X);
else
    gcamp = dff;
end
clear dff

%% Fit model for real and num_shuff iterations of shuffled data for each neuron


try
    %keyboard
    pc = parcluster('local');
    fprintf('%d',pc.NumWorkers)
    % set working area to local tmp file to avoid collisions
    %pc.JobStorageLocation = strcat('/tmp/',getenv('USER'),'-',getenv('SLURM_JOB_ID'));
    %start customized parallel pool
    pc.NumWorkers = 48;
    parpool(pc,48)
catch
    parpool(10);
end

switch statMethod
    case 'crossval'
        % Load or create data partition of trials (use same partitioning to allow for model comparison)
        try load(fullfile(fbasename,rec,sprintf('imagingReg_dataPartition_%dfolds_fs%d.mat',nCVFolds,frameRate)))
        catch
            load(fullfile(fbasename,rec,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']));
            
            ntrials = numel(tdtEvents.choice);
            trialStart = [1 tdtEvents.trialStart(2:end)-2*frameRate];
            trialStop  = [(tdtEvents.trialStart(2:end)-2*frameRate)-1 size(gcamp,1)];
            % Generate data partitions
            randGenerator = RandStream('mt19937ar','Seed',723176);
            
            dataPartition = cvpartition(ntrials,'KFold',nCVFolds,randGenerator);
            save(fullfile(fbasename,rec,sprintf('imagingReg_dataPartition_%dfolds_fs%d.mat',nCVFolds,frameRate)),'dataPartition','trialStart','trialStop')
        end
        
    case 'f'
        if reStart == 1
            load(fullfile(fbasename, rec, [fnext '.mat']), 'pvals','pvals_noRefit','F_shuff','F_nonShuff','F_shuff_noRefit','F_noRefit_nonShuff','con_iden','b','n','stats')
            roiStart = n+1;
            clear n
        else
            roiStart = 1;
        end

        for n = roiStart:size(gcamp,2)
            fprintf('\n%s%d\n', 'Calculating p-values for neuron ', n)
            % fit unshuffled data
            [F_nonShuff{n},F_noRefit_nonShuff{n},b{n},stats{n}]=get_fstat_eventReg(X,gcamp(:,n),preds_cells);
            % fit shuffled data
            temp_p = zeros(numel(preds_cells),num_shuff);
            temp_f = zeros(numel(preds_cells),num_shuff);
            parfor i=1:num_shuff
                gcamp_temp=circshift(gcamp(:,n),[randi(size(gcamp,1),1)]);
                gcamp_temp = gcamp_temp(:);
                [temp_f_refit(:,i),temp_f(:,i),temp_b(:,i)]=get_fstat_eventReg(X,gcamp_temp,preds_cells);
                fprintf('.')
            end
            F_shuff{n} = temp_f_refit;
            F_shuff_noRefit{n} = temp_f;
            b_shuff{n} = temp_b;
            
            % calculate p-values
            for g = 1:numel(preds_cells)
                pvals{n}(g) = sum(F_nonShuff{n}(g)<=F_shuff{n}(g,:))/num_shuff;
                pvals_noRefit{n}(g) = sum(F_noRefit_nonShuff{n}(g)<=F_shuff_noRefit{n}(g,:))/num_shuff;
            end
            save(fullfile(fbasename, rec, [fnext '.mat']), 'pvals','pvals_noRefit','F_shuff','F_nonShuff','F_shuff_noRefit','F_noRefit_nonShuff','con_iden','b','stats','n','-v7.3')            
        end
        delete(gcp('nocreate'))
end



