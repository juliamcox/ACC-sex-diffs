    function tdtEvents = checkTDTEvents(tdtEvents, data,allEvents) 


% find redundant events
temp = sort(allEvents(1,:));
reps = temp((diff(temp) == 0)); 
reps = unique(reps); 

% check each event type to make sure numbers match medpc data 

%trial start
if numel(tdtEvents.trialStart) ~= data.N
    if sum(data.raw.A ~=0) ~= data.N  %if session was force quit, probably an extra trial start time stamp (should also be in the medpc data), if so, exclude last timestamp
        tdtEvents.trialStart = tdtEvents.trialStart(1:end-1);
        tdtEvents.trialStartTimestamps = tdtEvents.trialStartTimestamps(1:end-1);
        tdtEvents.timestampError.trialStart = 1;
    end
    if numel(tdtEvents.trialStart)~=data.N
        % check if didn't stop recording
        idx = find(tdtEvents.trialStartTimestamps > tdtEvents.frameTimes_orig(end));
        if ~isempty(idx)
            tdtEvents.trialStart(idx) = [];
            tdtEvents.trialStartTimestamps(idx) = [];
        end
    end
    if numel(tdtEvents.trialStart)~=data.N
        % check for spurious event
        idxTrial = [];
        for i = 1:numel(reps)
            idxTrial = [idxTrial find(tdtEvents.trialStartTimestamps == reps(i))];
        end
        
        if isempty(idxTrial) || numel(idxTrial) > 1
            offFlag = 1;
            while offFlag
                if numel(tdtEvents.trialStart) ~= data.N
                    idxTrial = [];
                    if round(tdtEvents.trialStartTimestamps(1)-data.trialStart(1)) == round(tdtEvents.trialStartTimestamps(2)-data.trialStart(2))
                        offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
                    elseif  round(tdtEvents.nosePokeEntryTimestamps(1)-data.nosePokeEntry(1)) == round(tdtEvents.nosePokeEntryTimestamps(2)-data.nosePokeEntry(2))
                        offset = tdtEvents.nosePokeEntryTimestamps(1)-data.nosePokeEntry(1); %in case the first trial is the problem trial 
                    else
                        keyboard
                    end
                    temp = tdtEvents.trialStartTimestamps(1:numel(data.trialStart))-data.trialStart;
                    idxTrial = find(temp<offset,1,'first');
                    tdtEvents.trialStart(idxTrial) = [];
                    tdtEvents.trialStartTimestamps(idxTrial) = [];
                    tdtEvents.timestampError.trialStart = 1;
                else
                    offFlag = 0;
                end
            end
            
        else
            keyboard
            tdtEvents.trialStartTimestamps(idxTrial) = [];
            tdtEvents.trialStart(idxTrial) = [];
            tdtEvents.timestampError.trialStart = 1;
        end
    end
end



%nose poke entry, lever presentation
names = {'leverPresentation'; 'nosePokeEntry'};
for n = 1:numel(names)
    if numel(eval(sprintf('tdtEvents.%s', names{n}))) ~= data.N
        % check if didn't stop recording
        idx = find(eval(sprintf('tdtEvents.%sTimestamps',names{n})) > tdtEvents.frameTimes_orig(end));
        if ~isempty(idx)
           eval(sprintf('tdtEvents.%s(idx)=[];',names{n}));
           eval(sprintf('tdtEvents.%sTimestamps(idx)=[];',names{n}));
        end
    end
    if numel(eval(sprintf('tdtEvents.%s', names{n}))) ~= data.N
        % check for spurious event
        idxTrial = [];
        for i = 1:numel(reps)
            idxTrial = [idxTrial find(eval(sprintf('tdtEvents.%sTimestamps', names{n})) == reps(i))];
        end
        if numel(idxTrial) > 1
            idxTrial = [];
             offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
             repFlag = 1;
             while repFlag == 1
                 temp = eval(sprintf('tdtEvents.%sTimestamps(1:numel(data.%s))-data.%s;',names{n},names{n}, names{n}));
                 idxTrial = find(temp<offset,1,'first');
                 eval(sprintf('tdtEvents.%s(idxTrial) = [];', names{n}));
                 eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', names{n}));
                 if numel(eval(sprintf('tdtEvents.%s',names{n}))) == numel(eval(sprintf('data.%s',names{n})))
                     repFlag = 0;
                 end
             end
                             eval(sprintf('tdtEvents.timestampError.%s = 1;', names{n}));

        elseif isempty(idxTrial)
            keyboard
        else
            eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', names{n}));
            eval(sprintf('tdtEvents.%s(idxTrial) = [];', names{n}));
            eval(sprintf('tdtEvents.timestampError.%s = 1;', names{n}));
        end
    end
end


% lever press
if numel(tdtEvents.RLeverPress)~=sum(data.choice==0)
    % check if didn't stop recording
    idx = find(tdtEvents.RLeverPressTimestamps > tdtEvents.frameTimes_orig(end));
    if ~isempty(idx)
        tdtEvents.RLeverPress(idx) = [];
        tdtEvents.RLeverPressTimestamps(idx) = [];
    end
end
if numel(tdtEvents.RLeverPress) ~= sum(data.choice==0)
    
    
    if numel(tdtEvents.RLeverPress) > sum(data.choice==0)
        idxTrial = [];
        for i = 1:numel(reps)
            idxTrial = [idxTrial find(tdtEvents.RLeverPressTimestamps == reps(i))];
        end
        if isempty(idxTrial)
            keyboard
        elseif numel(idxTrial) > 1
            idxTrial = [];
            offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
            temp = tdtEvents.RLeverPressTimestamps(1:end-1)-data.leverPress(data.choice==0);
            idxTrial = find(temp<offset,1,'first');
              tdtEvents.RLeverPress(idxTrial) = [];
            tdtEvents.RLeverPressTimestamps(idxTrial) = [];
            tdtEvents.timestampError.RLever = 1; 
        else
            tdtEvents.RLeverPress(idxTrial) = [];
            tdtEvents.RLeverPressTimestamps(idxTrial) = [];
            tdtEvents.timestampError.RLever = 1; 
        end
    else
        keyboard
    end
end

if numel(tdtEvents.LLeverPress)~=sum(data.choice==1)
    % check if didn't stop recording
    idx = find(tdtEvents.LLeverPressTimestamps > tdtEvents.frameTimes_orig(end));
    if ~isempty(idx)
        tdtEvents.LLeverPress(idx) = [];
        tdtEvents.LLeverPressTimestamps(idx) = [];
    end
end

if numel(tdtEvents.LLeverPress) ~= sum(data.choice==1)
    if numel(tdtEvents.LLeverPress) > sum(data.choice==1)
        idxTrial = [];
         for i = 1:numel(reps)
            idxTrial = [idxTrial find(tdtEvents.LLeverPressTimestamps == reps(i))];
         end
         if isempty(idxTrial)
             keyboard
         elseif numel(idxTrial) > 1
             idxTrial = [];
             offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
             repFlag = 1;
             while repFlag == 1
                 temp = tdtEvents.LLeverPressTimestamps(1:sum(data.choice==1))-data.leverPress(data.choice==1);
                 idxTrial = find(temp<offset,1,'first');
                 eval(sprintf('tdtEvents.%s(idxTrial) = [];', 'LLeverPress'));
                 eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', 'LLeverPress'));
                 if numel(tdtEvents.LLeverPress) == numel(data.leverPress(data.choice==1))
                     repFlag = 0;
                 end
             end
             eval(sprintf('tdtEvents.timestampError.%s = 1;', 'LLever'));
         else
             tdtEvents.LLeverPress(idxTrial) = [];
             tdtEvents.LLeverPressTimestamps(idxTrial) = [];
             tdtEvents.timestampError.LLever = 1;
        end
    else
        keyboard
    end
end




% reward delivery, reward port entry, reward exit 
names = {'rewardDelivery'};
for n = 1:numel(names)
    if numel(tdtEvents.rewardDelivery)~=sum(data.reward==1)
        % check if didn't stop recording
        idx = find(tdtEvents.rewardDeliveryTimestamps > tdtEvents.frameTimes_orig(end));
        if ~isempty(idx)
            tdtEvents.rewardDelivery(idx) = [];
            tdtEvents.rewardDeliveryTimestamps(idx) = [];
        end
    end
    if eval(sprintf('numel(tdtEvents.%s) ~= sum(data.reward==1)', names{n}))
        if eval(sprintf('numel(tdtEvents.%s) > sum(data.reward==1)', names{n}))
            idxTrial = [];
            
            for i = 1:numel(reps)
                idxTrial = [idxTrial eval(sprintf('find(tdtEvents.%sTimestamps == reps(i))', names{n}))];
            end
            if isempty(idxTrial)
                keyboard
            elseif numel(idxTrial) > 1
                idxTrial = [];
                offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
                repFlag = 1;
                while repFlag == 1 
                    temp = tdtEvents.rewardDeliveryTimestamps(1:sum(data.reward==1))-data.rewardDelivery(data.reward==1);
                    idxTrial = find(temp<offset,1,'first');
                    eval(sprintf('tdtEvents.%s(idxTrial) = [];', names{n}));
                    eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', names{n}));
                    if numel(tdtEvents.rewardDeliveryTimestamps) == numel(data.rewardDelivery(data.reward==1))
                        repFlag = 0;
                    end                        
                end
                eval(sprintf('tdtEvents.timestampError.%s = 1;', names{n}));
            else
                eval(sprintf('tdtEvents.%s(idxTrial) = [];', names{n}));
                eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', names{n}));
                eval(sprintf('tdtEvents.timestampError.%s = 1;', names{n}));
            end
        else
            keyboard
        end
    end
end

% reward port entry and exit
if numel(tdtEvents.rewardEntry)~=sum(data.reward==1)
    % check if didn't stop recording
    idx = find(tdtEvents.rewardEntryTimestamps > tdtEvents.frameTimes_orig(end));
    if ~isempty(idx)
        tdtEvents.rewardEntry(idx) = [];
        tdtEvents.rewardEntryTimestamps(idx) = [];
    end
end
if numel(tdtEvents.rewardEntry) ~= sum(data.reward==1)
    if  numel(tdtEvents.rewardEntry) > sum(data.reward==1)
        idxTrial = [];
        for i = 1:numel(reps)
            idxTrial = [idxTrial find(tdtEvents.rewardEntryTimestamps == reps(i))];
        end
        if isempty(idxTrial)
            keyboard
        elseif numel(idxTrial) > 1
             idxTrial = [];
                offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
                repFlag = 1;
                while repFlag == 1 
                    temp = tdtEvents.rewardEntryTimestamps(1:sum(data.reward==1))-data.rewardEntry(data.reward==1);
                    idxTrial = find(temp<offset,1,'first');
                    eval(sprintf('tdtEvents.%s(idxTrial) = [];', 'rewardEntry'));
                    eval(sprintf('tdtEvents.%sTimestamps(idxTrial) = [];', 'rewardEntry'));
                    if numel(tdtEvents.rewardEntryTimestamps) == numel(data.rewardDelivery(data.reward==1))
                        repFlag = 0;
                    end                        
                end
                eval(sprintf('tdtEvents.timestampError.%s = 1;', 'rewardEntry'));
        else
            tdtEvents.rewardEntry(idxTrial) = [];
            tdtEvents.rewardEntryTimestamps(idxTrial) = [];
            tdtEvents.rewardExit(idxTrial) = [];
            tdtEvents.rewardExitTimestamps(idxTrial) = [];
            tdtEvents.timestampError.rewardEntry = 1;
        end
    elseif numel(tdtEvents.rewardEntry) ~= sum(data.rewardEntry ~= 0)
        tdtEvents.medRecoverFlag.rewardEntry = 1; 
    else %detection failure
        tdtEvents.trialExcludeReward = find(data.reward==1 & data.rewardEntry ==0);
        tdtEvents.trialExclude = tdtEvents.trialExcludeReward;
       % keyboard
    end
end

% CS plus and CS minus

if isfield(tdtEvents, 'CSPlus') && numel(tdtEvents.CSPlus) > 5
    if numel(tdtEvents.CSPlus)~=sum(data.reward==1)
        % check if didn't stop recording
        idx = find(tdtEvents.CSPlusTimestamps > tdtEvents.frameTimes_orig(end));
        if ~isempty(idx)
            tdtEvents.CSPlus(idx) = [];
            tdtEvents.CSPlusTimestamps(idx) = [];
        end
    end
    
    if numel(tdtEvents.CSPlus) ~= sum(data.reward==1)
        if numel(tdtEvents.CSPlus) > sum(data.reward==1)
            idxTrial = [];
            offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
            offFlag = 1;
            while offFlag
                if numel(tdtEvents.CSPlus) ~= sum(data.reward==1&data.choice~=-1)
                    temp = tdtEvents.CSPlusTimestamps(1:sum(data.reward==1&data.choice~=-1))-data.cue(data.reward==1&data.choice~=-1);
                    idxTrial = find(temp<offset,1,'first');
                    tdtEvents.CSPlus(idxTrial) = [];
                    tdtEvents.CSPlusTimestamps(idxTrial) = [];
                    tdtEvents.timestampError.CSPlus = 1;
                else
                    offFlag = 0;
                end
            end
            if isempty(idxTrial)
                keyboard
                %             elseif numel(idxTrial) > 1
                %                 keyboard
            
            end
        else
            keyboard
        end
    end
    if numel(tdtEvents.CSMinus)~=sum(data.reward==0)
        % check if didn't stop recording
        idx = find(tdtEvents.CSMinusTimestamps > tdtEvents.frameTimes_orig(end));
        if ~isempty(idx)
            tdtEvents.CSMinus(idx) = [];
            tdtEvents.CSMinusTimestamps(idx) = [];
        end
    end
    if numel(tdtEvents.CSMinus) ~= sum(data.reward==0&data.choice~=-1)
        if numel(tdtEvents.CSMinus) > sum(data.reward==0&data.choice~=-1)
            idxTrial = [];
            %             for i = 1:numel(reps)
            %                 idxTrial = [idxTrial find(tdtEvents.CSMinusTimestamps == reps(i))];
            %             end
            offset = tdtEvents.trialStartTimestamps(1)-data.trialStart(1);
            offFlag = 1;
            while offFlag
                if numel(tdtEvents.CSMinus) ~= sum(data.reward==0&data.choice~=-1)
                    temp = tdtEvents.CSMinusTimestamps(1:sum(data.reward==0&data.choice~=-1))-data.cue(data.reward==0&data.choice~=-1);
                    idxTrial = find(temp<offset,1,'first');
                    if isempty(idxTrial) 
                        if numel(tdtEvents.CSMinus) - sum(data.reward==0&data.choice~=-1) == 1
                            idxTrial = numel(tdtEvents.CSMinus)
                        else
                            keyboard
                        end
                    end
                            
                    tdtEvents.CSMinus(idxTrial) = [];
                    tdtEvents.CSMinusTimestamps(idxTrial) = [];
                    tdtEvents.timestampError.CSMinus = 1;
                else
                    offFlag = 0;
                end
            end
        else
            keyboard
        end
    end
    
    
end
