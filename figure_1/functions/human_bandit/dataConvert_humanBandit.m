function dataConvert_humanBandit(cohort)

% Data location
basefilename = fullfile(whereAreWe('figurecode'),'raw_data','figure_1','human_bandit',cohort); 
% Get data file names
cd(basefilename);
flist = dir('*.json'); 
flist = {flist(:).name};

for nf = 1:numel(flist)
    % if file has not been converted
  if ~exist([flist{nf}(1:end-5) '.mat'], 'file')
        % Read text
        data.raw = jsondecode(fileread(flist{nf}));
        % Find entries that contain trial epochs 
        idx = find(cellfun(@(x) isfield(x,'test_phase'),data.raw));
        
        % Find trial initiation latencies for all trials
        idx_ts                 = cellfun(@(x) contains(x.test_phase,'initiation'),data.raw(idx));
        idx_ts                 = idx(idx_ts); 
        idx_ts                 = idx_ts(5:end);
        % Find trials that were omitted at initiation if applicable (postpilot, trial timed out after 10 seconds of initiation screen)
        initOmitIdx            = cellfun(@(x) ~isfield(x,'rt'),data.raw(idx_ts)); 
        % initate data structure fill w/ NaN for omitted trials
        data.trialStartLatency       = zeros(size(idx_ts));
        data.leverJitter             = zeros(size(idx_ts));
        data.trialStart_cueOn        = zeros(size(idx_ts));
        data.trialStart_press        = zeros(size(idx_ts));
        
        data.trialStartLatency(initOmitIdx)= NaN;
        data.leverJitter(initOmitIdx)      = NaN;
        data.trialStart_cueOn(initOmitIdx) = NaN;
        data.trialStart_press(initOmitIdx) = NaN;

        
        idx_ts                 = idx_ts(~initOmitIdx); 
        data.initOmitIdx       = initOmitIdx; 
        
        trialStart_temp        = cellfun(@(x) x.rt, data.raw(idx_ts));
        leverJitter_temp       = cellfun(@(x) x.post_jitter, data.raw(idx_ts)); 
        trialStart_cueOn_temp  = cellfun(@(x) x.start_time, data.raw(idx_ts));
        trialStart_press_temp  = cellfun(@(x) x.key_time, data.raw(idx_ts)); 
        
        data.trialStartLatency(~initOmitIdx)  = trialStart_temp;
        data.leverJitter(~initOmitIdx) = leverJitter_temp;
        data.trialStart_cueOn(~initOmitIdx) = trialStart_cueOn_temp;
        data.trialStart_press(~initOmitIdx) = trialStart_press_temp;

      
        % Find choice, outcome and block 
        idx_ch                                    = cellfun(@(x) contains(x.test_phase,'choice'),data.raw(idx));
        idx_ch                                    = idx(idx_ch);
        idx_ch                                    = idx_ch(5:end);       
        idx_ch                                    = idx_ch(~initOmitIdx); %remove omitted initiation from index list
        
        tempChoice                                = cellfun(@(x) x.key_press, data.raw(idx_ch),'UniformOutput',false);
        data.choice                               = nan(size(initOmitIdx));
        data.reward                               = nan(size(initOmitIdx));
        data.blockID                              = nan(size(initOmitIdx));
        %Find omitted trials
        data.omitIdx                              = zeros(size(initOmitIdx));
        omitIdx                                   = cellfun(@isempty,tempChoice);
        data.omitIdx(~initOmitIdx)                = omitIdx; 
        data.omitIdx                              = logical(data.omitIdx); 
        
        % Find choice presentation
        data.choiceOn                             = nan(size(initOmitIdx)); 
        choiceOn                                  = cellfun(@(x) x.time_elapsed, data.raw(idx_ch)); 
        data.choiceOn(~initOmitIdx)               = choiceOn; 
        data.choiceOn_time                        = {nan(size(initOmitIdx))}; 
        choiceOn_time                             = cellfun(@(x) x.start_time, data.raw(idx_ch),'UniformOutput',false);
        data.choiceOn_time(~initOmitIdx)          = choiceOn_time; 
        data.choiceOn_time(data.omitIdx)          = {NaN};
        data.choiceOn_time(data.initOmitIdx)      = {NaN};
        data.choiceOn_time                        = cell2mat(data.choiceOn_time');
        
        % Find choice time
        data.choicePress                          = {nan(size(initOmitIdx))}; 
        choicePress                               = cellfun(@(x) x.key_time, data.raw(idx_ch),'UniformOutput',false);
        data.choicePress(~initOmitIdx)            = choicePress;
        data.choicePress(data.omitIdx)            = {NaN};
        data.choicePress(data.initOmitIdx)        = {NaN};
        data.choicePress                          = cell2mat(data.choicePress'); 
        
        % Find outcome presentation jitter
        data.outcomeJitter                        = nan(size(initOmitIdx)); 
        outcomeJitter                             = cellfun(@(x) x.post_jitter, data.raw(idx_ch)); 
        data.outcomeJitter(~initOmitIdx)          = outcomeJitter;
        % Choice latency
        data.choiceLatency                        = {nan(size(initOmitIdx))}; 
        choiceLatency                             = cellfun(@(x) x.rt, data.raw(idx_ch),'UniformOutput',false);
        data.choiceLatency(~initOmitIdx)          = choiceLatency;
        data.choiceLatency(data.omitIdx)          = {NaN};
        data.choiceLatency(data.initOmitIdx)      = {NaN};
        data.choiceLatency                        = cell2mat(data.choiceLatency');
        
        
        tempTempChoice                            = {nan(size(initOmitIdx))};
        tempTempChoice(~initOmitIdx)              = tempChoice;
        tempTempChoice(data.omitIdx)              = {NaN};
        tempTempChoice(data.initOmitIdx)          = {NaN};

        tempChoice                                = cell2mat(tempTempChoice);
         
        data.choice(tempChoice==65,1)             = 1;
        data.choice(tempChoice==76,1)             = 0;
        
        
        idx_ch                                    = cellfun(@(x) contains(x.test_phase,'choice'),data.raw(idx));
        idx_ch                                    = idx(idx_ch);
        idx_ch                                    = idx_ch(5:end);    
        data.reward                               = cellfun(@(x) double(x.rewarded), data.raw(idx_ch),'UniformOutput',false);
        data.reward(data.omitIdx)                 = {NaN};
        data.reward(data.initOmitIdx)             = {NaN};
        data.reward                               = cell2mat(data.reward); 
        tempBlock                                 = cellfun(@(x) x.correct_choice, data.raw(idx_ch), 'UniformOutput', false);
        data.blockID(contains(tempBlock,'left'),1)  = 1;
        data.blockID(contains(tempBlock,'right'),1) = 0;
        
        % Find outcome block 
        idx_out                                   = cellfun(@(x) contains(x.test_phase, 'outcome'),data.raw(idx));
        idx_out                                   = idx(idx_out);
        idx_out                                   = idx_out(5:end); 
        data.outcomeOn                            = cellfun(@(x) x.time_elapsed, data.raw(idx_out));
        data.outcomeOn(data.omitIdx)              = NaN;
        data.outcomeOn(data.initOmitIdx)          = NaN;
        
        save([flist{nf}(1:end-5) '.mat'], 'data')
        clear data
 end
end




