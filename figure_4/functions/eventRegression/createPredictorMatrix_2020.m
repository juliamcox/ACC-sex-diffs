function [x_basic,con_iden] = createPredictorMatrix_2020(rec,fnext,frameRate,frameRateO,regName,nframes,expt)

%%% Figure 2

% Generate predictor matrix for regression regName 

% This function has an incomplete section for continuous variables (i.e. position data) if add that in later

% rec:  location of the data
% fnext: file save name
% frameRate: desired frame rate for regression
% frameRateO: acquisition rate of the movie
% regName: which predictors? 
% nframes: length of recording
% expt: ACC_DMS_imaging or DMS imaging


%% Parameters  

% smoothing parameters for position data
smoothParams.method = 'median';
smoothParams.niter = 3;
smoothParams.winSize = 7;
smoothParams.nanHandle = 'omitnan';

% set pathnames
fbasename = fullfile(whereAreWe('imaging'),expt);
fbasename_bs = fullfile(whereAreWe('bucket'),'DMS_Bandit','basis sets'); 

% Get parameters for regName
[cons, con_shift, time_back_orig, time_forward_orig, type1, shift_con, bsIDs,contVar,fullVar] = getEvents(regName,frameRate);

%%%Extract behavioral events in samples
load(fullfile(fbasename, rec,'info.mat'), 'syncfilename', 'behavfilename'); %load filenames for behavior and sync files
 try load(fullfile(fbasename,rec,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']));
 catch
    fprintf('Extracting and downsampling behavior events \n')
    
    [tdtEvents] = behavEventRegression_2020(fullfile(fbasename, rec), syncfilename, behavfilename, frameRateO, 1000/round(frameRate));
    load(fullfile(fbasename,rec,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']));
 end

 %% Regressors 
 if ~isfield(tdtEvents,'ipsiChoice')
     try info = matfile(fullfile(fbasename, rec,'info.mat'));
         tdtEvents.ipsiChoice=tdtEvents.choice==info.ipsiSide;
         tdtEvents.ipsiChoice = double(tdtEvents.ipsiChoice);
         tdtEvents.ipsiChoice(tdtEvents.choice==-1) = -1;
         save(fullfile(fbasename,rec,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']),'tdtEvents','-append');
         
     catch
         keyboard
         save(fullfile(fbasename,rec,['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']),'tdtEvents','-append');
     end
 end
%Trial start
TrialStart = tdtEvents.trialStart;
%Nose poke
NPTimes = tdtEvents.nosePokeEntry;
NPTimesI = tdtEvents.nosePokeEntry(tdtEvents.ipsiChoice==1);
NPTimesC = tdtEvents.nosePokeEntry(tdtEvents.ipsiChoice==0);
NPTimesP = tdtEvents.nosePokeEntry(tdtEvents.reward==1);
NPTimesN = tdtEvents.nosePokeEntry(tdtEvents.reward==0&tdtEvents.choice~=-1);
idx = find(tdtEvents.reward==1)+1;
idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
NPPrevRew = tdtEvents.nosePokeEntry(idx);
idx = find(tdtEvents.reward==0&tdtEvents.choice~=-1)+1;
idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
NPPrevNRew = tdtEvents.nosePokeEntry(idx); 
%Lever presentation
LeverPresent = tdtEvents.leverPresentation;
LeverPresentI = tdtEvents.leverPresentation(tdtEvents.ipsiChoice==1);
LeverPresentC = tdtEvents.leverPresentation(tdtEvents.ipsiChoice==0);
LeverPresentP = tdtEvents.leverPresentation(tdtEvents.reward==1);
LeverPresentN = tdtEvents.leverPresentation(tdtEvents.reward==0&tdtEvents.choice~=-1);
%Lever press
LeverTimes = tdtEvents.leverPress;
LeverTimesI = tdtEvents.leverPress(tdtEvents.ipsiChoice==1);
LeverTimesC = tdtEvents.leverPress(tdtEvents.ipsiChoice==0);
LeverTimesP = tdtEvents.leverPress(tdtEvents.reward==1);
LeverTimesN = tdtEvents.leverPress(tdtEvents.reward==0&tdtEvents.choice~=-1);
%Outcome
CS = tdtEvents.outcome;
CSI = tdtEvents.outcome(tdtEvents.ipsiChoice==1);
CSC = tdtEvents.outcome(tdtEvents.ipsiChoice==0);
CS(tdtEvents.choice==-1) = []; %remove omitted trials
CSRew = tdtEvents.CSPlus;

rewIdx = find(tdtEvents.reward==1); 
rewIdx = setdiff(rewIdx,tdtEvents.trialExcludeReward); %in case of dropped frames/excluded trials
nrewIdx = find(tdtEvents.reward==0); 
nrewIdx = setdiff(nrewIdx,tdtEvents.trialExcludeNoReward); %in case of dropped frames/excluded trials
nrewIdx = nrewIdx(ismember(nrewIdx,find(tdtEvents.choice~=-1))); %remove omitted trials


CSRewI = tdtEvents.CSPlus(ismember(rewIdx,find(tdtEvents.ipsiChoice==1)));
CSRewC = tdtEvents.CSPlus(ismember(rewIdx,find(tdtEvents.ipsiChoice==0)));
CSNoRew = tdtEvents.CSMinus;
CSNoRewI = tdtEvents.CSMinus(ismember(nrewIdx,find(tdtEvents.ipsiChoice==1)));
CSNoRewC = tdtEvents.CSMinus(ismember(nrewIdx,find(tdtEvents.ipsiChoice==0)));
%Reward entry
RewardEnter = tdtEvents.rewardEntry;
RewardEnterI = tdtEvents.rewardEntry(ismember(rewIdx,find(tdtEvents.ipsiChoice==1)));
RewardEnterC = tdtEvents.rewardEntry(ismember(rewIdx,find(tdtEvents.ipsiChoice==0)));
%Reward exit
RewardExit = tdtEvents.rewardExit;
RewardExitI = tdtEvents.rewardExit(ismember(rewIdx,find(tdtEvents.ipsiChoice==1)));
RewardExitC = tdtEvents.rewardExit(ismember(rewIdx,find(tdtEvents.ipsiChoice==0))); 

% stay v switch (next trial) 
stay = [tdtEvents.choice(2:end)==tdtEvents.choice(1:end-1) NaN];
try
stay(find(tdtEvents.choice==-1)-1) = NaN;
catch
    omitIdx = find(tdtEvents.choice==-1);
    omitIdx = setdiff(omitIdx,1);
    stay(omitIdx-1) = NaN;
end
stay(find(tdtEvents.choice==-1)) = NaN;

CSRewStay = tdtEvents.CSPlus(ismember(rewIdx,find(stay==1)));
CSNoRewStay = tdtEvents.CSMinus(ismember(nrewIdx,find(stay==1))); 

%%% Make predictor matrix 

% Load basis set if necessary
if strcmp(type1,'spline')
    if iscell(bsIDs) 
        for nb = 1:numel(bsIDs)
            try load(fullfile(fbasename_bs, ['bs_' bsIDs{nb} '_' num2str(frameRate) '.mat']))
                bs{nb} = eval(bsIDs{nb}); 
            catch
                range = [1 1+str2double(bsIDs{nb}(1))*frameRate];
                idx = strfind(bsIDs{nb}, '_');
                numFun = str2double(bsIDs{nb}(idx(1)+1:idx(2)-1)); %number of functions in basis set
                order = 4; %order of bsplines (4 is cubic)
                basisobj = create_bspline_basis(range, numFun, order); %create basis function
                eval(sprintf('bs_%s =  getbasismatrix(range(1):range(2), basisobj);', bsIDs{nb})); %extract matrix of basis set
                save(fullfile(fbasename_bs, ['bs_' bsIDs{nb} '.mat']), ['bs_' bsIDs{nb}]);
            end
            basis_set{nb} = full(eval(sprintf('bs_%s',bsIDs{nb}))); 
        end
    else
        bs = bsIDs; % identify necessary basis sets
        try load(fullfile(fbasename_bs, ['bs_' bs '_' num2str(frameRate) '.mat']))
        catch
            range = [1 1+str2double(bs(1))*frameRate]; %number of time points that function will span
            idx = strfind(bs, '_');
            numFun = str2double(bs(idx(1)+1:idx(2)-1)); %number of functions in basis set
            order = 4; %order of bsplines (4 is cubic)
            basisobj = create_bspline_basis(range, numFun, order); %create basis function
            eval(sprintf('bs_%s =  getbasismatrix(1:1+str2num(bs(1))*frameRate, basisobj);', bs)); %extract matrix of basis set
            save(fullfile(fbasename_bs, ['bs_' bs '_' num2str(frameRate) '.mat']), ['bs_' bs]);
        end
        basis_set = full(eval(['bs_' bs]));
    end
end
x_basic = [];
con_iden = [];
event_times_mat = [];
% generate event vectors and convolve with basis sets if necessary 
for con = 1:numel(cons)
    if numel(time_back_orig)>1
        time_back = time_back_orig(con);
        time_forward = time_forward_orig(con); 
    else
        if con_shift(con) == 1 && shift_con == 1
            time_back = 0;
            time_forward = time_back_orig+time_forward_orig;
        else
            time_back = time_back_orig;
            time_forward = time_forward_orig;
        end
    end
   
   con_times = eval(cons{con});
   if isempty(strfind(cons{con}, 'I')) && isempty(strfind(cons{con}, 'C')) && isempty(strfind(cons{con}, 'Rew')) %all trial regressors
       con_times(tdtEvents.choice==-1) = []; %remove omitted trials
   end
   %con_times(isnan(con_times)) = []; % for choice regressors 
   % create vector where event times are 1 and other times 0
   con_binned = zeros(1,nframes);
   con_binned(int32(con_times))=1;
   
   if strcmp(type1,'spline') == 1
       con_binned = circshift(con_binned,-time_back*frameRate);
       event_times_mat = vertcat(event_times_mat,con_binned);
       if ~iscell(basis_set)
              for num_sets = 1:numel(basis_set(1,:))
               temp_conv_vec = conv(con_binned,basis_set(:,num_sets));
               x_basic = horzcat(x_basic,temp_conv_vec(1:numel(con_binned))');
           end
           con_iden = [con_iden ones(1,size(basis_set,2))*con];
       else
           for num_sets = 1:numel(basis_set{con}(1,:))
               temp_conv_vec = conv(con_binned,basis_set{con}(:,num_sets));
               x_basic = horzcat(x_basic,temp_conv_vec(1:numel(con_binned))');
           end
           con_iden = [con_iden ones(1,size(basis_set{con},2))*con];
       end
   elseif strcmp(type1,'time_shift')==1
        x_con=[];
        shift_back=frameRate*time_back;   %how many points do I shift forward and backwards
        shift_forward=frameRate*time_forward;

        for shifts = -shift_back:shift_forward
            x_con=horzcat(x_con,circshift(con_binned,[0,shifts])');
        end
        
        x_basic=horzcat(x_basic,x_con);
        con_iden=[con_iden ones(1,size(x_con,2))*con];
   end
end

if ~isempty(contVar)
    keyboard
    %%% extract tracking data 
    % load calibration data
    try load(fullfile(fbasename,rec, 'cal.mat'));
    catch
        keyboard
        cal = tdtTrackingCalibrations(rec);
        save(fullfile(fbasename,rec, 'cal.mat'),'cal');
    end
    % load TDT data
    load(fullfile(fbasename,  rec, syncfilename));
    % extract smoothed tracking data 
    tracking = extractTracking(data, smoothParams,cal,frameRate,motionTrack.inscopixStart,motionTrack.inscopixStop); 
    for ne = 1:numel(contVars)
        switch contVars{ne}
            case 'velocity'   
                x_basic = cat(2,x_basic,tracking.velocity);
            case 'acceleration'
                x_basic = cat(2,x_basic,tracking.acceleration);
            case 'position'
                x_basic = cat(2,x_basic,tracking.spatialBin);
        end        
    end
    % Full trial variables 
    for ne = 1:numel(fullVar)
        
        
    end
    
    
end


save(fullfile(fbasename,rec,['M_' fnext '.mat']), 'frameRate','bsIDs','x_basic','event_times_mat','con_iden','x_basic');