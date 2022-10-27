function extractData_opto(cohort,ext,qFile,zscoreFlag,binNum,savehere)
% Extract all opto data 

%% Parameters
valType = {'qChosenDiff';'qChosen';'qDiff';'qTot';'qLeft';'qRight'};


%% Animal IDs

% Generate male, female, opsin, yfp subject  lists 
aids_opsin_m = generateAnimalList([cohort '_male']);
if ~contains(cohort,'ACC')
    aids_yfp_m = generateAnimalList('DMS_yfp_male'); 
else
    aids_yfp_m = generateAnimalList([cohort '_yfp_male']);
end

ext_opsin_m = repmat({ext},size(aids_opsin_m));
ext_yfp_m = repmat({ext},size(aids_yfp_m));

aids_opsin_f = generateAnimalList([cohort '_female']);
if ~contains(cohort,'ACC')
    aids_yfp_f = generateAnimalList('DMS_yfp_female');
else
aids_yfp_f = generateAnimalList([cohort '_yfp_female']);
end

ext_opsin_f = repmat({ext},size(aids_opsin_f));
ext_yfp_f = repmat({ext},size(aids_yfp_f));

aids_all = cat(1,aids_opsin_m,aids_yfp_m,aids_opsin_f,aids_yfp_f);
%% Generate predictor matrix 

% Opsin: nphr 0, yfp 0
% Sex: male 0, female 1


Latency = [];
Sex = [];
Laser = [];
Opsin = [];
Value = cell(size(valType));
Animal = []; 
Trial = [];
Reward = []; 
PrevOutcome = [];
Choice = [];
Value_zscore = cell(size(valType));
Value_ptile = cell(size(valType));
Session = [];
SessionType = [];
Length = []; 
Stay = [];
High=[];
TrialStart = []; 
% Male, opsin
[animal,value,latency,laser,trials,prevoutcome,choice,value_zscore,value_ptile,session,sessionType,length,stay,high,trialstart] = extractData_group(aids_opsin_m,ext_opsin_m,valType,qFile,zscoreFlag,binNum);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,zeros(size(latency))'); 
Opsin = cat(1,Opsin,ones(size(latency))');
Laser = cat(1,Laser,laser');
Animal = cat(1,Animal,animal');
Trial = cat(1,Trial,trials'); 
PrevOutcome = cat(1,PrevOutcome,prevoutcome'); 
Choice = cat(1,Choice,choice'); 
High = cat(1,High,high');
TrialStart = cat(1,TrialStart,trialstart'); 
aCounter = max(animal); 
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
Session = cat(1,Session,session);
SessionType = cat(1,SessionType,sessionType');
Length = cat(1,Length,length');
Stay = cat(1,Stay,stay'); 
% Male, YFP
[animal,value,latency,laser,trials,prevoutcome,choice,value_zscore,value_ptile,session,sessionType,length,stay,high,trialstart] = extractData_group(aids_yfp_m,ext_yfp_m,valType,qFile,zscoreFlag,binNum);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,zeros(size(latency))'); 
Opsin = cat(1,Opsin,zeros(size(latency))');
Laser = cat(1,Laser,laser');
Animal = cat(1,Animal,(animal+aCounter)');
Trial = cat(1,Trial,trials'); 
PrevOutcome = cat(1,PrevOutcome,prevoutcome');
Choice = cat(1,Choice,choice'); 
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
aCounter = aCounter+max(animal); 
Session = cat(1,Session,session);
SessionType = cat(1,SessionType,sessionType');
Length = cat(1,Length,length');
Stay = cat(1,Stay,stay'); 
High = cat(1,High,high');
TrialStart = cat(1,TrialStart,trialstart');
% Female, opsin
[animal,value,latency,laser,trials,prevoutcome,choice,value_zscore,value_ptile,session,sessionType,length,stay,high,trialstart] = extractData_group(aids_opsin_f,ext_opsin_f,valType,qFile,zscoreFlag,binNum);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,ones(size(latency))'); 
Opsin = cat(1,Opsin,ones(size(latency))');
Laser = cat(1,Laser,laser');
Animal = cat(1,Animal,(animal+aCounter)');
Trial = cat(1,Trial,trials'); 
PrevOutcome = cat(1,PrevOutcome,prevoutcome'); 
Choice = cat(1,Choice,choice'); 
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
aCounter = aCounter+max(animal); 
Session = cat(1,Session,session);
SessionType = cat(1,SessionType,sessionType');
Length = cat(1,Length,length');
Stay = cat(1,Stay,stay'); 
High = cat(1,High,high');
TrialStart = cat(1,TrialStart,trialstart');

% Female, YFP
[animal,value,latency,laser,trials,prevoutcome,choice,value_zscore,value_ptile,session,sessionType,length,stay,high,trialstart] = extractData_group(aids_yfp_f,ext_yfp_f,valType,qFile,zscoreFlag,binNum);
% concatenate latency and predictors excluding session breaks (i.e. NaNs)
Latency = cat(1,Latency,latency');
Sex = cat(1,Sex,ones(size(latency))'); 
Opsin = cat(1,Opsin,zeros(size(latency))');
Laser = cat(1,Laser,laser');
Animal = cat(1,Animal,(animal+aCounter)');
Trial = cat(1,Trial,trials'); 
PrevOutcome = cat(1,PrevOutcome,prevoutcome'); 
Choice = cat(1,Choice,choice');
for nv = 1:numel(valType)
    Value{nv} = cat(1,Value{nv},value{nv}');
    Value_zscore{nv} = cat(1,Value_zscore{nv},value_zscore{nv}');
    Value_ptile{nv} = cat(1,Value_ptile{nv},value_ptile{nv}');
end
Session = cat(1,Session,session);
SessionType = cat(1,SessionType,sessionType');
Length = cat(1,Length,length');
Stay = cat(1,Stay,stay'); 
High = cat(1,High,high');
TrialStart = cat(1,TrialStart,trialstart');

% Remove NaNs (omitted trials)
idx = isnan(Value{1});

Latency(idx) = [];
Sex(idx) = [];
Opsin(idx) = [];
Laser(idx) = [];
Animal(idx) = [];
Trial(idx) = [];
PrevOutcome(idx) = []; 
Choice(idx) = []; 
Session(idx) = []; 
Length(idx) = [];
SessionType(idx) = [];
Stay(idx) = [];
High(idx) = [];
for nv = 1:numel(valType)
Value{nv}(idx) = [];
Value_zscore{nv}(idx) = [];
Value_ptile{nv}(idx) = [];
end


% remove session breaks
Latency(Laser==-10) = [];
Sex(Laser==-10) = [];
Opsin(Laser==-10) = [];
Animal(Laser==-10)=[];
Trial(Laser==-10) = [];
PrevOutcome(Laser==-10) = [];
Choice(Laser==-10) = [];
Session(Laser==-10) = [];
SessionType(Laser==-10) = [];
Length(Laser==-10) = [];
Stay(Laser==-10) = [];
High(Laser==-10)= [];
TrialStart(Laser==-10) = []; 

for nv = 1:numel(valType)
Value_zscore{nv}(Laser==-10) = [];
Value_ptile{nv}(Laser==-10) = [];
Value{nv}(Laser==-10)=[];
end
Laser(Laser==-10)=[];

Laser(Laser == 6) = 3;

if sum(Laser==4)~=0
    keyboard
    Latency(Laser==4) = [];
    Sex(Laser==4) = [];
    Opsin(Laser==4) = [];
    Animal(Laser==4)=[];
    Value(Laser==4)=[];
    Trial(Laser==4) = [];
    PrevOutcome(Laser==4) = [];
    Value_zscore(Laser==4) = [];
    Value_ptile(Laser==4) = [];
    Session(Laser==4) = [];
   
    Stay(Laser==4) = [];
    High(Laser==4) = [];
    Laser(Laser==4)=[];
    
    
end


Animal_all = Animal;
Latency_all = Latency;
Sex_all = Sex;
Opsin_all = Opsin;
Value_all = Value;
Trial_all = Trial;
PrevOutcome_all = PrevOutcome;
Choice_all = Choice; 
Value_zscore_all = Value_zscore;
Value_ptile_all = Value_ptile;
Session_all = Session;
Laser_all = Laser; 
Stay_all = Stay;
High_all = High;
TrialStart_all = TrialStart;
SessionType_all = SessionType; 
Length_all = Length; 


% Remove laser trials from control
% idx = find(Laser==2);
% idx = idx(Laser(idx-1)==0);
% Laser(idx-1)=-2;
% idx = find(Laser==1);
% idx = idx(Laser(idx-1)==0);
% Laser(idx-1)=-1; 
% idx = find(Laser==3);
% idx = idx(Laser(idx+1) == 0);
% Laser(idx+1) = -3;



Laser_outcome = zeros(size(Laser));
Laser_outcome(Laser==1) = 1;
Laser_outcome(Laser==2) = NaN;
Laser_outcome(Laser==3) = NaN;
Laser_np = zeros(size(Laser));
Laser_np(Laser==2) = 1;
Laser_np(Laser==1) = NaN;
Laser_np(Laser==3) = NaN; 
Laser_ITI = zeros(size(Laser));
Laser_ITI(Laser == 3)=1; 
Laser_ITI(Laser==1) = NaN;
Laser_ITI(Laser==2) = NaN;
%Laser_ITI(Laser==0) = 0;

% Laser_outcome = Laser==1;
% Laser_np = Laser==2;
% Laser_ITI = Laser==3;

save(fullfile(savehere, sprintf('%s_allData_opto_%s_zscore%d_bin%d_%s',cohort,ext,zscoreFlag,binNum,qFile)));