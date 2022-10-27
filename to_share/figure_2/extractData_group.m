function [aNum,val,latency,laser,trials,prevOutcome,choice,val_zscore,val_ptile,session,sessionType,length,stay,high,trialstart] = extractData_group(aids,ext,valType,qFile,zscoreFlag,binNum)
ptiles = [0 100/binNum:100/binNum:100];

basefilename = fullfile(whereAreWe('bucket'), 'Operant');
val = cell(size(valType));
val_zscore = cell(size(valType));
val_ptile = cell(size(valType));
aNum = [];
latency = [];
laser = [];
trials =[];
prevOutcome = [];
choice = [];
session = [];
sessionType = [];
length = [];
stay = [];
high = [];
trialstart = []; 
for na = 1:numel(aids)
    %% Extract laser sessions
    try
    load(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat'])); % Load behavior data
    catch
        keyboard
    end
    % If response times have not been calculated
    if ~isfield(data, 'rt')
        data = responseTimesGUI(data);
    end
    
    if zscoreFlag
        trialStart = data.rt.nosePoke_zscore;
    else
        trialStart = data.rt.nosePoke;
        trialStart(trialStart<0.01) = 0.01;
    end
    %
    %     %    Extract and concatenate value if necessary
    if ~isfield(data,'qFile')
        %      Load Q-values and extract laser sessions
        load(fullfile(basefilename, aids{na}, qFile));
        data = concatValue(data,qLearn,aids{na});
        data.qFile = qFile;
        save(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat']), 'data')
    elseif ~contains(data.qFile,qFile)
        % Load Q-values and extract laser sessions
        load(fullfile(basefilename, aids{na}, qFile));
        data = concatValue(data,qLearn,aids{na});
        data.qFile = qFile;
        save(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat']), 'data')
    end
    %
    
    
    
    % Calculate desired value
    for nv = 1:numel(valType)
        switch valType{nv}
            case 'qDiff'
                thisVal{nv} = data.QRight-data.QLeft;
            case 'qTot'
                thisVal{nv} = data.QRight+data.QLeft;
            case 'qChosen'
                thisVal{nv} = nan(size(data.choice));
                thisVal{nv}(data.choice==1) = data.QLeft(data.choice==1);
                thisVal{nv}(data.choice==0) = data.QRight(data.choice==0);
            case 'qChosenDiff'
                thisVal{nv} = nan(size(data.choice));
                thisVal{nv}(data.choice==1) = data.QLeft(data.choice==1)-data.QRight(data.choice==1);
                thisVal{nv}(data.choice==0) = data.QRight(data.choice==0)-data.QLeft(data.choice==0);
            case 'qRight'
                thisVal{nv} = data.QRight;
            case 'qLeft'
                thisVal{nv} = data.QLeft;
        end
    end
  
    
    for nv = 1:numel(valType)
        % Value zscored by session
        breaks = find(data.choice==-10);
        thisval_zscore{nv} = arrayfun(@(x,y) [nanzscore(thisVal{nv}(x:y)) -10 -10], [1 breaks(2:2:end-1)+1], [breaks(1:2:end)-1],'UniformOutput',false);
        thisval_zscore{nv} = cell2mat(thisval_zscore{nv});
    end
    
    % Find laser and non-laser trials
    %thisLaser = [0 data.laserType(1:end-1)];
    thisLaser = data.laserType;
    thisLaser(thisLaser==-1) = 0;
thisLaser(thisLaser==6) = 3;
laserIdx1          = find(thisLaser==1);
laserIdx2          = find(thisLaser==2);
laserIdx3          = find(thisLaser==3);
thisLaser(laserIdx1)= -1; % trial w/ laser
thisLaser(laserIdx2)= -2; % trial w/ laser
thisLaser(laserIdx2+1) = 2; % trial after laser
thisLaser(laserIdx1+1) = 1; % trial after laser
thisLaser(laserIdx3+1)= -3; % trial after laser (for ITI condition)
thisLaser(laserIdx3) = 3; 
%     thisLaser(thisLaser==-1) = 0;
%     ITI_idx =find(thisLaser==6);
%     thisLaser(ITI_idx) = 0;
%     thisLaser(ITI_idx-1) = 3; % shift ITI back to current trial
%     
    % trial counter
    breaks = find(data.choice==-10);
    thisTrial = arrayfun(@(x,y) [(x:y)-x+1 -10 -10], [1 breaks(2:2:end-1)+1],[breaks(1:2:end-1)-1],'UniformOutput',false);
    thisTrial = cell2mat(thisTrial);
    
    % trial start time 
    thisTrialStart = data.trialStart; 
    
    % previous outcome
    thisPrev = nan(size(thisTrial));
    thisPrev(find(data.reward==1)+1) = 1;
    thisPrev(find(data.reward==0)+1) = 0;
    thisStay = [NaN data.choice(1:end-1)==data.choice(2:end)];
    thisHigh = data.choice==(data.blocks.blockIDAll-1);
    idx = find(data.choice~=-10); % remove session borders
    
    % concatenate trials 
    %val  = cat(2,val,thisVal{nv}(idx));
    aNum = cat(2,aNum,ones(size(idx)).*na);
    latency = cat(2,latency,trialStart(idx));
    laser = cat(2,laser,thisLaser(idx));
    trials = cat(2,trials,thisTrial(idx));
    prevOutcome = cat(2,prevOutcome,thisPrev(idx));
    choice = cat(2,choice,data.choice(idx));
    stay = cat(2,stay,thisStay(idx));
    high = cat(2,high,thisHigh(idx));
    trialstart = cat(2,trialstart,thisTrialStart(idx)); 
    
    for nv = 1:numel(valType)
        val_zscore{nv} = cat(2,val_zscore{nv},thisval_zscore{nv}(idx));
    end
    sessionType = cat(2,sessionType,ones(size(idx)));
    length = cat(2,length,ones(size(idx)));
    
    % Extract date for each session and sort to create session predictor   
    flist = data.daylist;
    clear fdate
    for nf = 1:numel(flist)
        slashidx = strfind(flist{nf},'\');
        if isempty(slashidx)
            slashidx = strfind(flist{nf},'/');
        end
        try
        load(fullfile(whereAreWe('bucket'),'Operant',aids{na},flist{nf}(slashidx(end)+1:end)));
        catch
            load(fullfile(whereAreWe('bucket'),'Operant',aids{na},flist{nf}));
        end
            
        fdate(nf) = datenum(data.raw.header.Start_Date);
    end
    
    
    clear thisval_zscore  trialStart thisLaser thisTrial thisPrev thisVal{nv}_ptile thisSession
    
    if numel(trials)~=numel(laser)
        keyboard
    end
    
    for nv = 1:numel(valType)
        thisVal{nv} = thisVal{nv}(idx);
    end
    %% Extract non-laser sessions
    
    load(fullfile(basefilename, aids{na}, qFile));
    try
    fList = {qLearn.fList.name};
    catch
        fList = qLearn.fList;
    end
    fIdx = find(contains(fList,'.TAB'));
    fList = fList(contains(fList,'.TAB'));
    
    
    for nv = 1:numel(valType)
        switch valType{nv}
            case 'qDiff'
                Val_all{nv} = qLearn.QDiff;
            case 'qTot'
                Val_all{nv} = qLearn.QRight+qLearn.QLeft;
            case 'qChosen'
                Val_all{nv} = qLearn.QChosen;
            case 'qChosenDiff'
                Val_all{nv} = qLearn.QChosenDiff;
            case 'qRight'
                Val_all{nv} = qLearn.QRight;
            case 'qLeft'
                Val_all{nv} = qLearn.QLeft;
        end
    end
    for nf = 1:numel(fList)
        thisfile = sprintf('%s_TAB.mat',fList{nf}(2:end-4));
        % Load data file
        if exist(fullfile(basefilename,aids{na},thisfile)) ~= 0
            load(fullfile(whereAreWe('bucket'),'Operant',aids{na},thisfile));
            if round((datenum(data.raw.header.End_Time,'HH:MM:SS')-datenum(data.raw.header.Start_Time,'HH:MM:SS'))*24) == 2 || round((data.nosePokeEntry(end)-data.nosePokeEntry(1))/60)>90
%                  rew = nanmean(data.choice(data.idx.prevRew_all) == data.choice(data.idx.prevRew_all-1));
%                 nrew = nanmean(data.choice(data.idx.prevNRew_all) == data.choice(data.idx.prevNRew_all-1));
%                 
              %  if rew - nrew >= perfThresh && numel(data.choice)>200
               if numel(data.choice)>200
                thisLength = ones(size(data.choice));
                fdate = cat(2,fdate,datenum(data.raw.header.Start_Date));
                %if round((datenum(data.raw.header.End_Time,'HH:MM:SS')-datenum(data.raw.header.Start_Time,'HH:MM:SS'))*24) == 1 || round((data.nosePokeEntry(end)-data.nosePokeEntry(1))/60)<90
                %  thisLength = zeros(size(data.choice));
                %
                %                 end

                for nv = 1:numel(valType)
                    thisVal_temp{nv} = nan(size(data.choice));
                    thisVal_temp{nv} = Val_all{nv}(qLearn.session == fIdx(nf));
                    thisval_zscore{nv} = nanzscore(thisVal_temp{nv});
                end
                
                % Extract latencies
                if ~isfield(data, 'rt')
                    data = responseTimesGUI(data);
                end
                
                if zscoreFlag
                    trialStart = data.rt.nosePoke_zscore;
                else
                    trialStart = data.rt.nosePoke;
                    trialStart(trialStart==0) = 0.01;
                end
                
                % laser predictor (control sessions)
                thisLaser = zeros(size(data.choice));
                
                % Previous outcome
                thisPrev = nan(size(thisLaser));
                thisPrev(find(data.reward(1:end-1)==1)+1) = 1;
                thisPrev(find(data.reward(1:end-1)==0)+1) = 0;
                thisStay = [NaN data.choice(1:end-1)==data.choice(2:end)];
                thisHigh = data.choice==(data.blocks.blockIDAll-1);
                % Concatenate session's data
                aNum = cat(2,aNum,ones(size(thisLaser)).*na);
                latency = cat(2,latency,trialStart);
                laser = cat(2,laser,thisLaser);
                trials = cat(2,trials,1:numel(data.choice));
                prevOutcome = cat(2,prevOutcome,thisPrev);
                choice = cat(2,choice,data.choice);
                sessionType = cat(2,sessionType,zeros(size(data.choice)));
                stay = cat(2,stay,thisStay);
                high = cat(2,high,thisHigh);
                length = cat(2,length,thisLength);
                trialstart = cat(2,trialstart,data.trialStart);
                for nv = 1:numel(valType)
                    val_zscore{nv} = cat(2,val_zscore{nv},thisval_zscore{nv}');
                    thisVal{nv} = cat(2,thisVal{nv},thisVal_temp{nv}');
                end
                clear thisVal_temp trialStart thisLaser thisPrev data thisval_zscore thisVal_ptile thisLength thisStay
                end
            end
        else
            keyboard
        end
        
    end
%     
% %     % sort sessions by date to generate session predictor
    sorted = sort(fdate);
    sessionStart = find(aNum==na & trials==1); 
    sessionStop = [sessionStart(2:end)-1, find(aNum==na,1,'last')];
    thisSession = [];
    s=1;
    for ns = 1:numel(fdate)
        try
       thisSession = cat(1,thisSession,ones(size(sessionStart(ns):sessionStop(ns)))'.*(find(fdate(ns)==sorted))); 
        catch
            idx = (find(fdate(ns)==sorted));
            if s == 1
                thisSession = cat(1,thisSession,ones(size(sessionStart(ns):sessionStop(ns)))'.*idx(1));
                s = 2;
            elseif s == 2
                thisSession = cat(1,thisSession,ones(size(sessionStart(ns):sessionStop(ns)))'.*idx(2));
                s=1;
            end

        end
    end
%     
     session = cat(1,session,thisSession);
    try
    for nv = 1:numel(valType)
        tempBins = prctile(thisVal{nv}(trials(aNum==na)~=1&choice(aNum==na)~=-1),ptiles);
        [~,~,thisVal_ptile{nv}] = histcounts(thisVal{nv},tempBins);
        val_ptile{nv} = cat(2,val_ptile{nv},thisVal_ptile{nv});
        val{nv} = cat(2,val{nv},thisVal{nv});
    end
    catch
        keyboard
    end
end
end
