function valueExtraction_opto(aids, qFile,ext,perfThresh) 

% extracts and concatenates all control sessions (value, choice, latencies)
% that meet perfThresh (rew-no rew stay prob) and at least 200 trials for
% sessionLength ('long';'short';'both')

% excludes omitted trials 


basefilename = fullfile(whereAreWe('behavior'));

for na = 1:numel(aids)
    
    fprintf('Extracting %s data...\n', aids{na})
    
    % Find list of all sessions with extension ext
    allDir = dir(fullfile(basefilename, aids{na}, sprintf('*.%s',ext)));
    allDir = {allDir(:).name};
    
    % load value estimates
    load(fullfile(basefilename, aids{na}, qFile));
    % file list for q-values
    try
        thisFlist = {qLearn.fList.name};
    catch
        thisFlist = qLearn.fList;
    end
    
    flist = [];
    
    for nf = 1:numel(allDir)
        try
            
            load(fullfile(basefilename, aids{na}, [allDir{nf}(2:end-4) sprintf('_%s.mat',ext)])) % load data structure
            
            
            % extract trial type indices, if necessary
            if ~isfield(data, 'idx')
                data = getIndices(data);
                save(fullfile(basefilename, aids{na}, [allDir{nf}(2:end-4) sprintf('_%s.mat',ext)]),'data')
            end
            % If response times have not been calculated
            if ~isfield(data, 'rt')
                data = responseTimesGUI(data);
                save(fullfile(basefilename, aids{na}, [allDir{nf}(2:end-4) sprintf('_%s.mat',ext)]),'data')
            end
            %%% Calculate trialStart probability for each session following rewarded and unrewarded outcome
            rew = nanmean(data.choice(data.idx.prevRew_all) == data.choice(data.idx.prevRew_all-1));
            nrew = nanmean(data.choice(data.idx.prevNRew_all) == data.choice(data.idx.prevNRew_all-1));
            
            if numel(data.choice)>=200 && rew-nrew > perfThresh
                
                if nanmean(data.rt.nosePoke)>40
                    keyboard
                end
                
                fIdx = find(contains(thisFlist,allDir{nf}));
                if numel(fIdx)>1
                    keyboard
                    fIdx = fIdx(1);
                end
                if isempty(fIdx)
                    % If the session didn't get run with the Qlearning model
                    warning(sprintf('No Q-values for session %s', allDir{nf}));
                    
                else
                    
                    
                    trialIdx = find(qLearn.session == fIdx); % indices for this session's q-values
                    if numel(trialIdx) ~= numel(data.choice)
                        keyboard
%                         % Extract response times for trials following non-omitted trials (omitted trials excluded for fitting q-learning model)
%                         %%% find omitted trials
%                         omitIdx   = find(data.choice==-1);
%                         
%                         if ~isempty(omitIdx)
%                             
%                             thisRight     = nan(size(data.choice))';
%                             thisLeft      = nan(size(data.choice))';
%                             
%                             thisChosen    = nan(size(data.choice))';
%                             thisUnchosen  = nan(size(data.choice))';
%                             
%                             qRight    = qLearn.QRight(trialIdx);
%                             qLeft     = qLearn.QLeft(trialIdx);
%                             qChosen(qLearn.choice(trialIdx)==1,1) = qLeft(qLearn.choice(trialIdx)==1);
%                             qChosen(qLearn.choice(trialIdx)==0,1) = qRight(qLearn.choice(trialIdx)==0);
%                             qUnchosen(qLearn.choice(trialIdx)==0,1) = qLeft(qLearn.choice(trialIdx)==0);
%                             qUnchosen(qLearn.choice(trialIdx)==1,1) = qRight(qLearn.choice(trialIdx)==1);
%                             
%                             trialCount = 1;
%                             trialCount2 = 1;
%                             try
%                                 for no = 1:numel(omitIdx)
%                                     
%                                     if omitIdx(no)+1 == numel(data.choice) || omitIdx(no) == numel(data.choice)
%                                         thisIdx = trialCount:omitIdx(no)-1;
%                                         thisRight(thisIdx) = qRight(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisLeft(thisIdx) = qLeft(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisChosen(thisIdx) = qChosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisUnchosen(thisIdx) = qUnchosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         trialCount = omitIdx(no)+2;
%                                         trialCount2 = trialCount2+numel(thisIdx);
%                                     else
%                                         thisIdx = trialCount:omitIdx(no);
%                                         thisRight(thisIdx) = qRight(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisLeft(thisIdx) = qLeft(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisChosen(thisIdx) = qChosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         thisUnchosen(thisIdx) = qUnchosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                                         trialCount = omitIdx(no)+2;
%                                         trialCount2 = trialCount2+numel(thisIdx);
%                                     end
%                                 end
%                             catch
%                                 keyboard
%                             end
%                             thisIdx = trialCount:numel(thisRight);
%                             thisRight(thisIdx) = qRight(trialCount2:trialCount2+numel(thisIdx)-1);
%                             thisLeft(thisIdx) = qLeft(trialCount2:trialCount2+numel(thisIdx)-1);
%                             thisChosen(thisIdx) = qChosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                             thisUnchosen(thisIdx) = qUnchosen(trialCount2:trialCount2+numel(thisIdx)-1);
%                             
%                             thisIdx = find(data.choice~=-1);
%                             thisIdx = thisIdx(ismember(thisIdx,find(~isnan(thisLeft)))); % remove trial after the omitted trial as well as omitted trial
%                             QRight_temp{nf} = thisRight(thisIdx);
%                             QLeft_temp{nf} = thisLeft(thisIdx);
%                             QChosen_temp{nf} = thisChosen(thisIdx);
%                             QUnchosen_temp{nf} = thisUnchosen(thisIdx);
%                             
%                             trialStart_temp_zscore{nf} = data.rt.nosePoke_zscore(thisIdx);
%                             leverPress_temp_zscore{nf} = data.rt.leverPress_zscore(thisIdx);
%                             withdraw_temp_zscore{nf}   = data.rt.withdraw_zscore(thisIdx);
%                             
%                             trialStart_temp{nf} = data.rt.nosePoke(thisIdx);
%                             leverPress_temp{nf} = data.rt.leverPress(thisIdx);
%                             withdraw_temp{nf}   = data.rt.withdraw(thisIdx);
%                             
%                             choice_temp{nf} = data.choice(thisIdx);
%                             ipsiChoice_temp{nf}  = data.ipsiChoice(thisIdx);
%                             stay_temp{nf} = data.choice(thisIdx)==data.choice(thisIdx-1);
%                             high_temp{nf}        = data.choice(thisIdx)==data.blocks.blockIDAll(thisIdx)-1;
%                             prevOutcome_temp{nf} = data.reward(thisIdx-1);
%                             session_temp{nf}     = ones(size(thisIdx)).*fIdx;
%                             laserType_temp{nf}        = data.laserType(thisIdx); 
%                         else %%% extract value for this session
%                             thisRight = qLearn.QRight(trialIdx);
%                             thisLeft  = qLearn.QLeft(trialIdx);
%                             qChosen(qLearn.choice(trialIdx)==1,1) = thisLeft(qLearn.choice(trialIdx)==1);
%                             qChosen(qLearn.choice(trialIdx)==0,1) = thisRight(qLearn.choice(trialIdx)==0);
%                             qUnchosen(qLearn.choice(trialIdx)==0,1) = thisLeft(qLearn.choice(trialIdx)==0);
%                             qUnchosen(qLearn.choice(trialIdx)==1,1) = thisRight(qLearn.choice(trialIdx)==1);
%                             
%                             thisIdx = find(data.choice~=-1);
%                             QRight_temp{nf} = thisRight(thisIdx);
%                             QLeft_temp{nf} = thisLeft(thisIdx);
%                             QChosen_temp{nf} = qChosen(thisIdx);
%                             QUnchosen_temp{nf} = qUnchosen(thisIdx);
%                             
%                             
%                             trialStart_temp_zscore{nf} = data.rt.nosePoke_zscore(thisIdx);
%                             leverPress_temp_zscore{nf} = data.rt.leverPress_zscore(thisIdx);
%                             withdraw_temp_zscore{nf}   = data.rt.withdraw_zscore(thisIdx);
%                             
%                             trialStart_temp{nf}  = data.rt.nosePoke(thisIdx);
%                             leverPress_temp{nf}  = data.rt.leverPress(thisIdx);
%                             withdraw_temp{nf}    = data.rt.withdraw(thisIdx);
%                             
%                             choice_temp{nf}      = data.choice(thisIdx);
%                             ipsiChoice_temp{nf}  = data.ipsiChoice(thisIdx);
%                             stay_temp{nf}        = data.choice(thisIdx)==data.choice(thisIdx-1);
%                             high_temp{nf}        = data.choice(thisIdx)==data.blocks.blockIDAll(thisIdx)-1;
%                             prevOutcome_temp{nf} = data.reward(thisIdx-1);
%                             session_temp{nf}     = ones(size(thisIdx)).*fIdx;
%                             laserType_temp{nf}        = data.laserType(thisIdx);
%                         end
%                         
                        
                    else
                        %%% extract value for this session
                        thisRight = qLearn.QRight(trialIdx);
                        thisLeft  = qLearn.QLeft(trialIdx);
                        qChosen = nan(size(thisRight));
                        qUnchosen = nan(size(thisRight));
                        qChosen(qLearn.choice(trialIdx)==1,1) = thisLeft(qLearn.choice(trialIdx)==1);
                        qChosen(qLearn.choice(trialIdx)==0,1) = thisRight(qLearn.choice(trialIdx)==0);
                        qUnchosen(qLearn.choice(trialIdx)==0,1) = thisLeft(qLearn.choice(trialIdx)==0);
                        qUnchosen(qLearn.choice(trialIdx)==1,1) = thisRight(qLearn.choice(trialIdx)==1);
                        
                       % thisIdx = find(data.choice~=-1);
                        thisIdx = 1:numel(data.choice);
                        
                        QRight_temp{nf} = [thisRight(thisIdx); -10];
                        QLeft_temp{nf} = [thisLeft(thisIdx); -10];
                        QChosen_temp{nf} = [qChosen(thisIdx); -10];
                        QUnchosen_temp{nf} = [qUnchosen(thisIdx); -10];
                        
                        
                        trialStart_temp_zscore{nf} = [data.rt.nosePoke_zscore(thisIdx) -10];
                        leverPress_temp_zscore{nf} = [data.rt.leverPress_zscore(thisIdx) -10];
                        withdraw_temp_zscore{nf}   = [data.rt.withdraw_zscore(thisIdx) -10];
                        
                        trialStart_temp{nf}  = [data.rt.nosePoke(thisIdx) -10];
                        leverPress_temp{nf}  = [data.rt.leverPress(thisIdx) -10];
                        withdraw_temp{nf}    = [data.rt.withdraw(thisIdx) -10];
                        
                        choice_temp{nf}      = [data.choice(thisIdx) -10];
                        ipsiChoice_temp{nf}  = [data.ipsiChoice(thisIdx) -10];
                        try
                            stay_temp{nf}        = [double(data.choice(thisIdx)==data.choice(thisIdx-1)) -10];
                            prevOutcome_temp{nf} = [double(data.reward(thisIdx-1)) -10];
                        catch
                            stay_temp{nf}        = [NaN data.choice(thisIdx(2:end))==data.choice(thisIdx(2:end)-1) -10];
                            prevOutcome_temp{nf} = [NaN data.reward(thisIdx(2:end)-1) -10];
                        end
                        high_temp{nf}        = [data.choice(thisIdx)==data.blocks.blockIDAll(thisIdx)-1 -10];
                       
                        session_temp{nf}     = [ones(size(thisIdx)).*fIdx -10];
                        laserType_temp{nf}    = [data.laserType(thisIdx) -10];
                        trialTime_temp{nf}   = [data.trialStart(thisIdx) -10];

                    end
                    
                    % ipsi/contra
                    if contains(data.laserSide, 'L')
                        QIpsi_temp{nf} = [QLeft_temp{nf}; -10];
                        QContra_temp{nf} = [QRight_temp{nf}; -10];
                    elseif contains(data.laserSide, 'R')
                        QIpsi_temp{nf} = [QRight_temp{nf}; -10];
                        QContra_temp{nf} = [QLeft_temp{nf}; -10];
                    else
                        keyboard
                    end
                    
                    
                    try
                        flist = cat(1,flist,allDir(nf));
                    catch
                        keyboard
                    end
                end
            end
            omit{nf} = sum(data.choice==-1);
            totalTrial{nf} = numel(data.choice);
            
            
        catch
            keyboard
        end
    end
    if exist('QLeft_temp') == 0
        qLeft =[];
        qRight = [];
        qChosen = [];
        qUnchosen = [];
        qIpsi  = [];
        qContra = [];
        qDiff = [];% relative value
        qTot =  []; % total value
        qChosenDiff = [];
        
        % Latencies
        trialStart = [];
        withdraw = [];
        leverPress = [];
        
        trialStart_zscore = [];
        withdraw_zscore = [];
        leverPress_zscore = [];
        
        stay = [];
        leftChoice = [];
        prevOutcome = [];
        idx = [];
        session = [];
        highProb = [];
        stay = [];
        leftChoice = [];
        prevOutcome = [];
        session = [];
        omit = [];
        totalTrial = [];
        laserType = [];
    else
        qLeft = cell2mat(QLeft_temp');
        qRight = cell2mat(QRight_temp');
        qChosen = cell2mat(QChosen_temp');
        qUnchosen = cell2mat(QUnchosen_temp');
        qIpsi  = cell2mat(QIpsi_temp');
        qContra = cell2mat(QContra_temp');
        qDiff = qRight-qLeft; % relative value
        qTot =  qRight+qLeft; % total value
        qChosenDiff = qChosen-qUnchosen;
        
        % Latencies
        trialStart = cell2mat(trialStart_temp);
        withdraw = cell2mat(withdraw_temp);
        leverPress = cell2mat(leverPress_temp);
        
        trialStart_zscore = cell2mat(trialStart_temp_zscore);
        withdraw_zscore = cell2mat(withdraw_temp_zscore);
        leverPress_zscore = cell2mat(leverPress_temp_zscore);
        
        idx = cellfun(@(x) ~isempty(x), stay_temp);     
        stay = cell2mat(stay_temp(idx));
        highProb = cell2mat(high_temp(idx));
        leftChoice = cell2mat(choice_temp);
        prevOutcome = cell2mat(prevOutcome_temp);
        session = cell2mat(session_temp); 
        omit = cell2mat(omit);
        totalTrial = cell2mat(totalTrial); 
        laserType = cell2mat(laserType_temp);
        trialTime = cell2mat(trialTime_temp);
    end
    if ~contains(qFile, '.mat')
        qFile = cat(2,qFile,'.mat');
    end
    save(fullfile(basefilename,aids{na},sprintf('value%s_perfThresh%s_%s',ext,num2str(perfThresh),qFile)),'trialTime','laserType','highProb','omit','totalTrial','qLeft','qRight',...
        'qChosen','qUnchosen','qIpsi','qContra','qDiff','qTot','qChosenDiff','trialStart','withdraw','leverPress','trialStart_zscore','withdraw_zscore','leverPress_zscore','idx','stay','leftChoice','prevOutcome','flist','session')
    save(fullfile(basefilename,aids{na},sprintf('value%s_flist_perfThresh%s_%s',ext,num2str(perfThresh),qFile)),'flist')

end
end
