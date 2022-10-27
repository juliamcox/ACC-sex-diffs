function allCtrlData(aids_f,aids_m,sessionLength,perfThresh,qFile,savehere)

% Parameters

valType = {'qChosenDiff';'qChosen';'qDiff';'qTot';'qLeft';'qRight'};


%% Extract data from all subjects
Latency = [];
Sex = [];
Value = cell(size(valType));
Animal = [];
Trial = [];
PrevOutcome = [];
Choice = [];
Value_zscore = cell(size(valType));
Value_ptile = cell(size(valType));
Session = [];
Stay = [];
% Male, opsin
[animal,value,latency,trials,prevoutcome,choice,value_zscore,value_ptile,session,stay] = extractData_ctrl(aids_m,valType,qFile,sessionLength,perfThresh);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,zeros(size(latency))');
Animal = cat(1,Animal,animal');
Trial = cat(1,Trial,trials');
PrevOutcome = cat(1,PrevOutcome,prevoutcome');
Choice = cat(1,Choice,choice');
aCounter = max(animal);
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
Session = cat(1,Session,session');
Stay = cat(1,Stay,stay');


% Female
[animal,value,latency,trials,prevoutcome,choice,value_zscore,value_ptile,session,stay] = extractData_ctrl(aids_f,valType,qFile,sessionLength,perfThresh);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,ones(size(latency))');
Animal = cat(1,Animal,(animal+aCounter)');
Trial = cat(1,Trial,trials');
PrevOutcome = cat(1,PrevOutcome,prevoutcome');
Choice = cat(1,Choice,choice');
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
Session = cat(1,Session,session');
Stay = cat(1,Stay,stay');


% Remove NaNs (omitted trials)
idx = isnan(Value{1});

Latency(idx) = [];
Sex(idx) = [];
Animal(idx) = [];
Trial(idx) = [];
PrevOutcome(idx) = [];
Choice(idx) = [];
Session(isnan(Session)) = [];
Stay(idx) = [];
for nv = 1:numel(valType)
    Value{nv}(idx) = [];
    Value_zscore{nv}(idx) = [];
    Value_ptile{nv}(idx) = [];
end


Latency(Latency<0.01) = 0.01;
save(fullfile(savehere,sprintf('ctrlSesssions_%s_perfThresh_%d_%s',sessionLength,perfThresh,qFile)),'Latency','Sex','Animal','Trial','PrevOutcome','Choice','Session','Stay','Value','Value_zscore','Value_ptile','aids_m','aids_f');
