function  [tdtEvents] = behavEventRegression_2020(filebasename, syncfilename, behavfilename, inscopixFs, binsize)

%Extracts event and motion tracking data from TDT, extracts lick times from
%medpc code (only if timestamps not messed up by lasers running in another box pre 4/2018), timestamps given in seconds (inscopix clock) and frames 

load(fullfile(filebasename, syncfilename)); % load TDT data
info = matfile(fullfile(filebasename, 'info.mat')); 
try 
    
    info.Properties.Writable = true;
    info.nframesLog = info.nFramesLog;
    info.exposureTime = 1/info.frameRate;
catch
end
   

tdtEvents.trialExcludeReward = [];
tdtEvents.trialExcludeNoReward = [];
tdtEvents.trialExclude = [];
medFlag = 0;


%Make list of timestamps in seconds
timestamps = linspace(0,length(data.streams.gcmp.data)/data.streams.gcmp.fs, length(data.streams.gcmp.data));
%Find inscopix frame times from TDT
tdtEvents.frameIdx = [];
inscopixOn = (find(data.streams.gcmp.data > 3.5)); 
z = diff(inscopixOn); 
tdtEvents.frameIdx(1) = inscopixOn(1);
tdtEvents.frameIdx = [tdtEvents.frameIdx inscopixOn(find(z>1)+1)]; % Index of inscopix frame onset at data.streams.gcmp.fs
if length(tdtEvents.frameIdx) ~= info.nframesLog +info.droppedCount
    tdtEvents.frameIdx = [];
    inscopixOn = (find(data.streams.gcmp.data > 4));
    z = diff(inscopixOn);
    tdtEvents.frameIdx(1) = inscopixOn(1);
    tdtEvents.frameIdx = [tdtEvents.frameIdx inscopixOn(find(z>1)+1)]; % Index of inscopix frame onset at data.streams.gcmp.fs
%     if length(tdtEvents.frameIdx)-1 == info.nframesLog + info.droppedCount %if there's one extra timestamp it is most likely a spurious end frame 
%         tdtEvents.frameIdx = tdtEvents.frameIdx(1:end-1); 
    if length(tdtEvents.frameIdx) ~= info.nframesLog +info.droppedCount
        warning('Number of frames wrong in behavEventRegression_2020')
        keyboard
    end
else
    keyboard %manual check of timestamps w/ new inscopix bug
end

tdtTime = (timestamps(tdtEvents.frameIdx(end))+(info.exposureTime/1000)-timestamps(tdtEvents.frameIdx(1)));

% compare TDT timestamps to inscopix computer clock 
recDur = info.recDur;
ratio = double(recDur/vpa(tdtTime(1)));

if (ratio > 1.01 || ratio < .98)
    warning('Timing mismatch')
   % keyboard
    try load(fullfile(filebasename, 'behavSync.mat'), 'ratio')
    catch ME
        fprintf(ME.message)
    end
end

%timestamps = (timestamps.*ratio); 

tdtEvents.frameTimes = timestamps(tdtEvents.frameIdx); % Times of inscopix frame onset in seconds
tdtEvents.tdtFs = data.streams.gcmp.fs; %acquisition frame rate for the frame times
tdtEvents.frameIdx_tdt = tdtEvents.frameIdx;
tdtEvents.frameIdx = []; 




% downsample frame times (will not do anything if already in desired sampling frequency)
endPoint = (size(tdtEvents.frameTimes,2)-1)*(1000/inscopixFs);
currTaxis = linspace(0,endPoint,size(tdtEvents.frameTimes,2))';
newTaxis = (0:binsize:endPoint)';
if newTaxis(end)<endPoint; newTaxis = [newTaxis; endPoint]; end
tdtEvents.frameTimes_orig = tdtEvents.frameTimes;
tdtEvents.frameTimes = interp1q(currTaxis,tdtEvents.frameTimes',newTaxis);
z = diff(tdtEvents.frameTimes); 
if abs(1000/binsize - (1/mean(z))) > .1
    warning('Frame detection error')
    keyboard
end

%TDT events
event_names = {'rewardDelivery'; 'rewardEntry'; 'RLeverPress'; 'LLeverPress'; 'trialStart'; 'nosePokeEntry'; 'leverPresentation'; 'CSPlus'; 'CSMinus'};
tdtEvents.event_names = [event_names; {'rewardExit';'allNose'; 'lickTimes'; 'nosePokeExit'}];

% Extract event times and frames for event_names
allEventsTDT = [];
try
for i = 1:length(event_names)
    temp = getfield(data.epocs, ['up0' num2str(i)]);
    eventTimes = temp.onset'; %Time of the events in seconds
    eventFrames = [];
    for ii = 1:length(temp.onset)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - temp.onset(ii))); 
    end
    tdtEvents = setfield(tdtEvents,[event_names{i} 'Timestamps'], eventTimes);
    tdtEvents = setfield(tdtEvents,[event_names{i}], eventFrames);
    allEventsTDT = [allEventsTDT [eventTimes; ones(size(eventTimes)).*i]];
    clear temp eventTimes eventTimestamps eventFrames
end
catch ME
    medFlag = 1;
end

% Find reward exit time
tdtEvents.rewardExitTimestamps = data.epocs.dn02.onset';
for ii = 1:length(data.epocs.dn02.onset)
    [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - data.epocs.dn02.onset(ii))); 
end
tdtEvents.rewardExit = eventFrames; 
allEventsTDT = [allEventsTDT [tdtEvents.rewardExitTimestamps; ones(size(tdtEvents.rewardExitTimestamps)).*(i+1)]]; 

%if recorded all nose poke entries
if isfield(data.epocs,'up10') 
    eventFrames = [];
    eventTimes = data.epocs.up10.onset';
    for ii = 1:length(eventTimes)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - eventTimes(ii)));
    end
    tdtEvents.allNoseTimestamps = eventTimes;
    tdtEvents.allNoseTimes = eventFrames; 
    % Check whether lick times were ever up10????
    if mean(eventTimes-data.epocs.dn10.onset')~=-0.01
        keyboard
    end
    
end

%if recorded lick times with tdt (after repaired cable)
if isfield(data.epocs,'up13') 
    eventFrames = [];
    eventTimes = data.epocs.up13.onset';
    for ii = 1:length(eventTimes)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - eventTimes(ii)));
    end
    tdtEvents.lickTimestamps = eventTimes;
    tdtEvents.lickTimes = eventFrames; 
end

%if recorded motif
if isfield(data.epocs,'up14') 
    eventFrames = [];
    eventTimes = data.epocs.up14.onset';
    for ii = 1:length(eventTimes)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - eventTimes(ii)));
    end
    tdtEvents.cameraTimestamps = eventTimes;
    tdtEvents.cameraTimes = eventFrames; 
end


%% Timestamps from medpc

% MedPC events for double check and lick times 
load(fullfile(filebasename, behavfilename));


% check TDT for spurious events 
tdtEvents = checkTDTEvents(tdtEvents, data,allEventsTDT); 

% use first and last trial start to align med times
tdtTime = tdtEvents.trialStartTimestamps(1,end)-tdtEvents.trialStartTimestamps(1,1);
temp = data.trialStart(data.trialStart~=0); %in case of mis-writing of medpc data
medTime = temp(end)-data.trialStart(1); 

ratio2 = double(vpa(tdtTime)/vpa(medTime));

if ratio2 > 1.01 || ratio2 < 0
    warning('Timing mismatch')
   % keyboard
end


firstEvent = tdtEvents.trialStartTimestamps(1); %seconds
firstEventMedPC = data.trialStart(1); 
offset = firstEvent-(firstEventMedPC*ratio2);

if ~isfield(data, 'rightCue')
    data.rightCue = data.cue;
end

% lick times (if not recorded in tdt, convert from medpc)
clear eventTimestamps

% generate a time-varying offset vector to convert lick times (to deal with occasional
% slowing down of medpc clock in old recordings
tempOffset1 = tdtEvents.trialStartTimestamps(1:data.N) - (data.trialStart(1:data.N).*ratio2);
tempOffset2 = tdtEvents.nosePokeEntryTimestamps - data.nosePokeEntry.*ratio2;
temp = sort(([tdtEvents.RLeverPressTimestamps tdtEvents.LLeverPressTimestamps])) - data.leverPress(data.choice~=-1).*ratio2;
tempOffset3(data.choice==-1) = NaN;
tempOffset3(data.choice~=-1) = temp;
try %some recordings the CS sync signal did not work
    tempOffset4 = sort([tdtEvents.CSPlusTimestamps tdtEvents.CSMinusTimestamps]) - data.rightCue.*ratio2;
catch
    tempOffset4 = [];
end

tempOffset = nanmean([tempOffset1; tempOffset2; tempOffset3; tempOffset4]);
tempOffsetLick = [tempOffset1 tempOffset2 tempOffset3 tempOffset4]; %combine offset estimates from all events for lick times
tempOffsetLickAxis = [data.trialStart(1:data.N) data.nosePokeEntry(1:data.N) data.leverPress(1:data.N)]; %times of all the events for lick offset estimate
if ~isempty(tempOffset4)
    tempOffsetLickAxis = [tempOffsetLickAxis data.rightCue];
end
[~,sortIdx] = sort(tempOffsetLickAxis,'ascend');
tempOffsetLickAxis = tempOffsetLickAxis(sortIdx)';
tempOffsetLick = tempOffsetLick(sortIdx)';


[minTS,lickTS1] = min(abs(tempOffsetLickAxis - data.raw.Z(1))); %find closest event to first time

lickTSend = find(data.raw.Z <= tempOffsetLickAxis(end),1,'last');
lickOffset = interp1q(tempOffsetLickAxis(lickTS1:end), tempOffsetLick(lickTS1:end), linspace(data.raw.Z(1),data.raw.Z(lickTSend),numel(data.raw.Z(1:lickTSend)))');
lickOffset = [lickOffset; repmat(lickOffset(end), numel(data.raw.Z)-numel(lickOffset),1)];

idx = find(isnan(lickOffset)); % happens if the first lick timestamp is before the first trial start

for n = flipud(idx)'
    lickOffset(n) = lickOffset(n+1);
end



if isfield(tdtEvents,'lickTimes')
    tdtEvents.lickTimestamps_medPC = (data.raw.Z.*ratio2)' + lickOffset;
    eventFrames = [];
    clear eventFrames
    for ii = 1:numel(tdtEvents.lickTimestamps_medPC)
        if isnan(tdtEvents.lickTimestamps_medPC(ii))
            eventFrames(ii) = NaN;
        else
            [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - tdtEvents.lickTimestamps_medPC(ii)));
        end
    end
    tdtEvents.lickTimes_medPC = eventFrames;
else
    tdtEvents.lickTimestamps = (data.raw.Z.*ratio2)' + lickOffset;
    eventFrames = [];
    clear eventFrames
    for ii = 1:numel(tdtEvents.lickTimestamps)
        if isnan(tdtEvents.lickTimestamps(ii))
            eventFrames(ii) = NaN;
        else
            [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - tdtEvents.lickTimestamps(ii)));
        end
    end
    tdtEvents.lickTimes = eventFrames;
    
end

% detect lick bouts that are separated by at least 1 second 
d = diff(tdtEvents.lickTimestamps); %in seconds
offsetIdx = find(d>=1); 
onsetIdx = find(d>=1)+1; 
c = tdtEvents.lickTimes(offsetIdx(2:end))-tdtEvents.lickTimes(onsetIdx(1:end-1));
try
onsetIdx = [1; onsetIdx];
offsetIdx = [offsetIdx; numel(tdtEvents.lickTimestamps)]; 

catch
    
onsetIdx = [1; onsetIdx'];
offsetIdx = [offsetIdx'; numel(tdtEvents.lickTimestamps)]; 
end
tdtEvents.lickOffset = tdtEvents.lickTimes(offsetIdx);
tdtEvents.lickOnset = tdtEvents.lickTimes(onsetIdx); 

% convert medPC nose poke exit to tdt timestamps and add time varying
% offset
clear eventFrames
tdtEvents.nosePokeExitTimestamps = data.nosePokeExit.*ratio2 + tempOffset; 
for ii = 1:length(tdtEvents.nosePokeExitTimestamps)
    [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - tdtEvents.nosePokeExitTimestamps(ii)));
end
tdtEvents.nosePokeExit = eventFrames;

% medpc events
event_names = {'trialStart'; 'nosePokeEntry'; 'nosePokeExit'; 'leverPresentation'; 'leverPress'; 'rewardDelivery'; 'rewardEntry'; 'rewardExit'};
medEvents.event_names = [event_names; {'cueTimes'}];
allEvents = []; 

for i = 1:length(event_names)
    eventTimes = getfield(data, event_names{i}); % - firstEventMedPC;
    eventTimes = eventTimes(eventTimes ~= 0); 
    eventTimes = eventTimes.*ratio2 + offset; 
    for ii = 1:length(eventTimes)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - eventTimes(ii))); %Frame index of events 
    end
    medEvents = setfield(medEvents,[event_names{i} 'Timestamps'], eventTimes);
    medEvents = setfield(medEvents,[event_names{i}], eventFrames);
    allEvents = [allEvents; eventTimes' ones(size(eventTimes')).*i];
    if sum(eventTimes<0) > 0 
        fprintf('\n\nNegative timestamp: %s', event_names{i})
        %keyboard
    end
    clear temp eventTimes eventFrames
end


if medFlag
    % TDT didn't acquire cue times, check accuracy of medTimes and if offset too large, perform trial by trial alignment
    keyboard
    tdtEvents.CSPlusTimestamps = data.rightCue(data.reward==1).*ratio2 + tempOffset(data.reward==1);
    tdtEvents.CSMinusTimestamps = data.rightCue(data.reward==0).*ratio2 + tempOffset(data.reward==0);
    clear eventFrames
    for ii = 1:numel(tdtEvents.CSPlusTimestamps)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - tdtEvents.CSPlusTimestamps(ii)));
    end
    tdtEvents.CSPlus = eventFrames;
    clear eventFrames
    for ii = 1:numel(tdtEvents.CSMinusTimestamps)
        [~,eventFrames(ii)] = min(abs(tdtEvents.frameTimes - tdtEvents.CSMinusTimestamps(ii)));
    end
    tdtEvents.CSMinus = eventFrames;
    %Test
    %         testNosePoke = (data.nosePokeEntry.*ratio2) + tempOffset;
    %         testLeverPress = (data.leverPress.*ratio2) + tempOffset;
    %         testcue = (data.rewardDelivery(data.reward==1).*ratio2) + tempOffset(data.reward==1);
    %
    %         figure(), plot(tdtEvents.nosePokeEntryTimestamps-testNosePoke);
    %         figure(), plot(sort([tdtEvents.RLeverPressTimestamps tdtEvents.LLeverPressTimestamps]) - testLeverPress);
    %         figure(), plot(tdtEvents.rewardDeliveryTimestamps - testcue);
    
end


%dropped frames

if info.droppedCount ~= 0
    keyboard %update from orig version
    warning('There are dropped frames in this recording')
    droppedIdx = info.droppedIdx;
   % droppedIdx = droppedIdx.*((1000/binsize)/inscopixFs);
    x = 1;
    xx = 1;
    for i = 1:length(data.reward)
        if data.reward(i) == 1
            tdtEvents.trialTimeReward(x,:) = [tdtEvents.trialStart(i)-2*(1000/binsize), tdtEvents.rewardExit(x)+3*(1000/binsize), i];
            x=x+1;
        elseif data.reward(i) == 0 & data.choice(i)~=-1
            tdtEvents.trialTimeNoReward(xx,:) = [tdtEvents.trialStart(i)-2*(1000/binsize), tdtEvents.CSMinus(xx)+.5+4*(1000/binsize), i];
            xx = xx+1;
        end
    end
    
    %find trials that contain dropped frames
    
    for i = 1:length(info.droppedIdx)-1
       % keyboard
        if droppedIdx(i+1)-droppedIdx(i) ~= 1
            [temp,~] = find(droppedIdx(i) >= tdtEvents.trialTimeReward(:,1) & droppedIdx(i) < tdtEvents.trialTimeReward(:,1));
            if ~isempty(temp)
                tdtEvents.trialExcludeReward = [tdtEvents.trialExcludeReward temp];
            end
            clear temp
            [temp,~] = find(droppedIdx(i) >= tdtEvents.trialTimeNoReward(:,1) & droppedIdx(i) < tdtEvents.trialTimeNoReward(:,1));
            if ~isempty(temp)
                tdtEvents.trialExcludeNoReward = [tdtEvents.trialExcludeNoReward temp];
            end
            clear temp
        end
    end
    i = i+1;
    [temp,~] = find(droppedIdx(i) >= tdtEvents.trialTimeReward(:,1) & droppedIdx(i) <= tdtEvents.trialTimeReward(:,2));
    tdtEvents.trialExcludeReward = [tdtEvents.trialExcludeReward temp];
    clear temp
    [temp,~] = find(droppedIdx(i) >= tdtEvents.trialTimeNoReward(:,1) & droppedIdx(i) <= tdtEvents.trialTimeNoReward(:,2));
    tdtEvents.trialExcludeNoReward = [tdtEvents.trialExcludeNoReward temp];
    clear temp
    
    
    tdtEvents.trialExcludeNoReward = unique(tdtEvents.trialExcludeNoReward);
    tdtEvents.trialExcludeReward = unique(tdtEvents.trialExcludeReward);
    tdtEvents.trialExclude = [tdtEvents.trialExclude tdtEvents.trialTimeReward(tdtEvents.trialExcludeReward,3)' tdtEvents.trialTimeNoReward(tdtEvents.trialExcludeNoReward,3)'];

end


tdtEvents.reward = data.reward;
tdtEvents.choice = data.choice;
try
    tdtEvents.ipsiChoice = data.choice==info.ipsiSide;
    tdtEvents.ipsiChoice = double(tdtEvents.ipsiChoice); 
    tdtEvents.ipsiChoice(tdtEvents.choice==-1) = -1; 
catch
    keyboard
end
 
tdtEvents.rewardIdx = find(data.reward == 1); 
tdtEvents.rewardIdx2 = find(data.rewardEntry>0); 
if length(tdtEvents.rewardIdx) ~= length(tdtEvents.rewardIdx2)
    warning('Missed reward port entries')
end
tdtEvents.leverPress = NaN.*zeros(size(tdtEvents.choice));
tdtEvents.leverPress(tdtEvents.choice==0) = tdtEvents.RLeverPress(~isnan(tdtEvents.RLeverPress));
tdtEvents.leverPress(tdtEvents.choice==1) = tdtEvents.LLeverPress(~isnan(tdtEvents.LLeverPress));


%outcome
temp = sort([(tdtEvents.CSMinus) (tdtEvents.rewardDelivery)]);
tdtEvents.outcome(data.choice~=-1) = temp;
tdtEvents.outcome(data.choice==-1) = NaN; 

save(fullfile(filebasename,['behavEventsRegression_bs' num2str(binsize) '.mat']), 'tdtEvents', 'ratio', 'ratio2')



