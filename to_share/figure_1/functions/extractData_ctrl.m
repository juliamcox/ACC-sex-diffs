function [aNum,val,latency,trials,prevOutcome,choice,val_zscore,val_ptile,session,stay] = extractData_ctrl(aids,valType,qFile,sessionLength,perfThresh)
binNum =4;
ptiles = [0 100/binNum:100/binNum:100];
basefilename = fullfile(whereAreWe('bucket'), 'Operant');
val = cell(size(valType));
val_zscore = cell(size(valType));
val_ptile = cell(size(valType));
aNum = [];
latency = [];
trials =[];
prevOutcome = []; 
choice = [];
session = [];

stay = [];
for na = 1:numel(aids)
    
    %% Extract laser sessions
    % Load value and behavior data 
    try qLearn=load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)));
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist; 
    catch
        valueExtraction_TAB(aids(na),qFile,sessionLength,perfThresh);
        qLearn=load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)));
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist;
    end
    
   
    
    % Calculate desired value
    for nv = 1:numel(valType)
        thisVal{nv} = eval(sprintf('qLearn.%s',valType{nv}))';
    end
    
    sessionList = unique(qLearn.session);

    for nv = 1:numel(valType)
        % Value zscored by session
        thisval_zscore{nv} = arrayfun(@(x) nanzscore(thisVal{nv}(qLearn.session==x)),sessionList,'UniformOutput',false);
        thisval_zscore{nv} = cell2mat(thisval_zscore{nv});
    end
    
    
    % trial counter 
    thisTrial = arrayfun(@(x) 1:sum(qLearn.session==x), sessionList,'UniformOutput',false);
    thisTrial = cell2mat(thisTrial);
    
    
    % previous outcome
    thisPrev = qLearn.prevOutcome;
    
    %val  = cat(2,val,thisVal{nv}(idx));
    aNum = cat(2,aNum,ones(size(thisPrev)).*na); 
    latency = cat(2,latency,qLearn.trialStart);
    trials = cat(2,trials,thisTrial);
    prevOutcome = cat(2,prevOutcome,thisPrev);
    choice = cat(2,choice,qLearn.leftChoice); 
    stay = cat(2,stay,qLearn.stay); 
    for nv = 1:numel(valType)
        val_zscore{nv} = cat(2,val_zscore{nv},thisval_zscore{nv});
    end

    session = cat(2,session,qLearn.session); 
    
    for nv = 1:numel(valType)
        tempBins = prctile(thisVal{nv},ptiles);
        [~,~,thisVal_ptile{nv}] = histcounts(thisVal{nv},tempBins);
        val_ptile{nv} = cat(2,val_ptile{nv},thisVal_ptile{nv});
        val{nv} = cat(2,val{nv},thisVal{nv});
    end
    
    
    clear thisval_zscore  trialStart  thisTrial thisPrev thisVal_ptile thisSession

    
end
end
