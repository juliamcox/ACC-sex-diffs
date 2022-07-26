function extractData(ids_m,ids_f, perfThresh,qFile,weightSpreadsheet,ext_laser,binNum,binNum_choice,cohort,sessionLength,thresh)
savehere = fullfile(whereAreWe('figurecode'),'processed_data'); 
[~,sheetnames]=xlsfinfo(weightSpreadsheet); %sheet names of dominance spreadsheet which has all of the weights
behaviorTable = [];
for na = 1:numel(ids_m)
    fprintf('Processing %s \n', ids_m{na})
    [thisBehavior, flist_m{na}] = extractAllBehavior(ids_m{na},perfThresh,qFile,sheetnames,weightSpreadsheet,ext_laser,0,binNum,binNum_choice,cohort,sessionLength,thresh);
    
    behaviorTable = cat(1,behaviorTable,thisBehavior);
end
for na = 1:numel(ids_f)
    fprintf('Processing %s \n', ids_f{na})
    
    [thisBehavior, flist_f{na}] = extractAllBehavior(ids_f{na},perfThresh,qFile,sheetnames,weightSpreadsheet,ext_laser,1,binNum,binNum_choice,cohort,sessionLength,thresh);
    
    behaviorTable = cat(1,behaviorTable,thisBehavior);
end

% add in opsin
behaviorTable.opsin = zeros(size(behaviorTable.laser));
aids = cat(1,generateAnimalList(sprintf('%s_female',cohort)),generateAnimalList(sprintf('%s_male',cohort)));
idx = cell2mat(cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids,'UniformOutput',false));
behaviorTable.opsin(idx) = 1; 

behaviorTable.trialInit(behaviorTable.trialInit<0.01,:) = 0.01; % standardize min latency w/ different acquisition frequencies
behaviorTable.weight_zscore = nanzscore(behaviorTable.weight); 

% find average weight per subject
behaviorTable.meanWeight = zeros(size(behaviorTable,1),1);
aids = cat(1,ids_m,ids_f);
for na = 1:numel(aids)
   idx = find(strcmp(behaviorTable.aID, aids{na}));
   sessionIdx = unique(behaviorTable.session(idx));
   sessionIdx(isnan(sessionIdx)) = [];
   meanWeight = arrayfun(@(x) behaviorTable.weight(find(behaviorTable.session==x&strcmp(behaviorTable.aID, aids{na}),1,'first')),sessionIdx);
   behaviorTable.meanWeight(idx) = nanmean(meanWeight); 
end
behaviorTable.meanWeight_zscore = nanzscore(behaviorTable.meanWeight); 

save(fullfile(savehere,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort,sessionLength, num2str(thresh),qFile)),'behaviorTable','flist_f','flist_m')