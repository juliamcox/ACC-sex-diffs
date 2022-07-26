function extractQParams(aids_f,aids_m, perfThresh, sessionLength, qFile,ext,savehere)



params_f = extract(aids_f,qFile,perfThresh,sessionLength,ext);
params_m = extract(aids_m, qFile,perfThresh,sessionLength,ext);

save(fullfile(savehere,sprintf('qParams_perfThres%d_%s_%s_%s.mat', perfThresh,sessionLength,qFile,ext)),'params_f','params_m')

end

function params = extract(aids, qFile, perfThresh, sessionLength, ext)

fbasename = fullfile(whereAreWe('bucket'),'Operant'); % data location

for na = 1:numel(aids)
    cd(fullfile(fbasename, aids{na}));
    fprintf('Processing %s...\n', aids{na})
    
    % load value estimates
    load(fullfile(fbasename, aids{na}, qFile));
    
    % initialize variables for trialStart probability
    thisAlpha  = [];
    thisBeta  = [];
    thisStay  = [];
    thisSide  = [];
    if contains(qFile,'2alpha')
        thisAlpha2 = [];
    end
    % Find list of all sessions with extension ext
    allDir = dir(['*.', ext]);
    allDir = {allDir(:).name};
    
    
    sessionCounter = 1;
    for nf = 1:numel(allDir)
        try load([allDir{nf}(2:end-4) '_' ext '.mat']) % load data structure
        
        
        % Check session length
        if (contains(sessionLength,'short') && (round((datenum(data.raw.header.End_Time,'HH:MM:SS')-datenum(data.raw.header.Start_Time,'HH:MM:SS'))*24) == 1 || round((data.nosePokeEntry(end)-data.nosePokeEntry(1))/60)<90)) || (contains(sessionLength,'long') && (round((datenum(data.raw.header.End_Time,'HH:MM:SS')-datenum(data.raw.header.Start_Time,'HH:MM:SS'))*24) == 2 || round((data.nosePokeEntry(end)-data.nosePokeEntry(1))/60)>90)) || contains(sessionLength,'both')
            
            % extract trial type indices, if necessary
            if ~isfield(data, 'idx')
                data = getIndices(data);
                save([allDir{nf}(2:end-5) '.mat'])
            end
            
            %%% Calculate trialStart probability for each session following rewarded and unrewarded outcome
            rew = nanmean(data.choice(data.idx.prevRew_all) == data.choice(data.idx.prevRew_all-1));
            nrew = nanmean(data.choice(data.idx.prevNRew_all) == data.choice(data.idx.prevNRew_all-1));
            
            if rew - nrew >= perfThresh
                
                %%% Generate predictor matrix with value and trial number,animalID and session
                try
                    valIdx = find(contains({qLearn.fList(:).name}, allDir{nf}));
                catch
                    valIdx = find(contains(qLearn.fList, allDir{nf}));
                end
                
                if ~isempty(valIdx)
                    thisAlpha = cat(1,thisAlpha,qLearn.alpha_session(valIdx));
                    thisBeta  = cat(1,thisBeta,qLearn.beta_session(valIdx));
                    thisStay  = cat(1, thisStay,qLearn.stay_session(valIdx));
                    thisSide  = cat(1, thisSide, qLearn.side_session(valIdx));
                    if contains(qFile,'2alpha')
                        thisAlpha2 = cat(1,thisAlpha2,qLearn.alpha2_session(valIdx));
                    end
                end  
            end
        end
        catch
            keyboard
        end
    end
    
    params.alpha(na)   = mean(thisAlpha);
    params.beta(na)    = mean(thisBeta);
    params.stay(na)    = mean(thisStay);
    params.side(na)    = mean(thisSide);
    params.alpha_sem(na)   = std(thisAlpha)./sqrt(numel(thisAlpha));
    params.beta_sem(na)    = std(thisBeta)./sqrt(numel(thisBeta));
    params.stay_sem(na)    = std(thisStay)./sqrt(numel(thisBeta));
    params.side_sem(na)    = std(thisSide)./sqrt(numel(thisSide));
    params.alpha_count(na)   = (numel(thisAlpha));
    params.beta_count(na)    = (numel(thisBeta));
    params.stay_count(na)    = (numel(thisBeta));
    params.side_count(na)    = (numel(thisSide));
    if contains(qFile,'2alpha')
        params.alpha2(na)  = mean(thisAlpha);
    end
    try
        params.alpha_m(na) = qLearn.alpha;
        params.beta_m(na)  = qLearn.beta;
        params.stay_m(na)  = qLearn.stay;
        params.side_m(na)  = qLearn.bias;
        if contains(qFile,'2alpha')
            params.alpha2_m(na)  = qLearn.alpha2;
        end
    catch
        params.alpha_m(na) = qLearn.alpha_mouse;
        params.beta_m(na)  = qLearn.beta_mouse;
        params.stay_m(na)  = qLearn.stay_mouse;
        params.side_m(na)  = qLearn.bias_mouse;
         if contains(qFile,'2alpha')
            params.alpha2_m(na)  = qLearn.alpha2_mouse;
        end
    end
    
    
    
end
end
