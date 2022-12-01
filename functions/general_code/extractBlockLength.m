function extractBlockLength

%% Short sessions
% female subjects
aids = generateAnimalList('ACC_DMS_nphr_female');
aids = cat(1,aids,generateAnimalList('ACC_DMS_nphr_yfp_female'));
aids = cat(1,aids,generateAnimalList('DMS_nphr_d1_female'));
aids = cat(1,aids,generateAnimalList('DMS_nphr_d2_female'));
aids = cat(1,aids,generateAnimalList('DMS_yfp_female')); 
aids = cat(1,aids, generateAnimalList('imaging_female'));
% male subjects
aids = cat(1,aids,generateAnimalList('ACC_DMS_nphr_male'));
aids = cat(1,aids,generateAnimalList('ACC_DMS_nphr_yfp_male'));
aids = cat(1,aids,generateAnimalList('DMS_nphr_d1_male'));
aids = cat(1,aids,generateAnimalList('DMS_nphr_d2_male'));
aids = cat(1,aids,generateAnimalList('DMS_yfp_male')); 
aids = cat(1,aids, generateAnimalList('imaging_male'));

sessionLength = 'short';
qFile         = 'qLearn_session_all_2022.mat';
perfThresh    = 0.1;

[blockNum,blockNum_avg,blockNum_sem,blockNum_std,blockLen]=extractSession(aids,perfThresh,qFile,sessionLength);
[blockNum_e,blockNum_avg_e,blockNum_sem_e,blockNum_std_e,blockLen_e]=extractSession(generateAnimalList('estrous'),perfThresh,'qLearn_session_estrous_2022.mat',sessionLength);

blockNum = cat(2,blockNum,blockNum_e);
blockNum_avg = cat(2,blockNum_avg,blockNum_avg_e);
blockNum_sem = cat(2,blockNum_sem,blockNum_sem_e);
blockNum_std = cat(2,blockNum_std,blockNum_std_e);
blockLen = cat(2,blockLen,blockLen_e);

sessionLength = 'long';

[blockNum_long,blockNum_avg_long,blockNum_sem_long,blockNum_std_long,blockLen_long]=extractSession(aids,perfThresh,qFile,sessionLength);

[blockNum_e_long,blockNum_avg_e_long,blockNum_sem_e_long,blockNum_std_e_long,blockLen_e_long]=extractSession(generateAnimalList('estrous'),perfThresh,'qLearn_session_estrous_2022.mat',sessionLength);


blockNum_long = cat(2,blockNum_long,blockNum_e_long);
blockNum_avg_long = cat(2,blockNum_avg_long,blockNum_avg_e_long);
blockNum_sem_long = cat(2,blockNum_sem_long,blockNum_sem_e_long);
blockNum_std_long = cat(2,blockNum_std_long,blockNum_std_e_long);
blockLen_long = cat(2,blockLen_long,blockLen_e_long);

keyboard

save(fullfile(whereAreWe('figurecode'),'processed_data','blockLength.mat'), 'blockNum','blockNum_avg','blockNum_sem','blockNum_std','blockLen','blockNum_long','blockNum_avg_long','blockNum_sem_long','blockNum_std_long','blockLen_long'); 


keyboard
end


function [blockNum,blockNum_avg,blockNum_sem,blockNum_std,blockLen_all,blockLen]=extractSession(aids,perfThresh,qFile,sessionLength)
basefilename = fullfile(whereAreWe('behavior'));
blockLen = cell(size(aids));
 blockLen_all = [];
for na = 1:numel(aids)
    
    fprintf('Processing %s...\n', aids{na})
    
    try
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist;
    catch
        valueExtraction_TAB(aids(na),qFile,sessionLength,perfThresh);
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist;
    end
   
    for ns = 1:numel(flist)
        try
       load(fullfile(basefilename,aids{na}, [flist{ns}(2:end-4),'_TAB.mat']))
        
       try
       blockNum{na}(ns) = numel(data.blocks.trialIdx);
       blockLen{na} = cat(2,blockLen{na},data.blocks.lenDist(1:end-1));
       blockLen_all = cat(2,blockLen_all,data.blocks.lenDist(1:end-1));
       catch
       blockNum{na}(ns) = numel(data.blocks.switch_idx);
       blockLen{na} = cat(2,blockLen{na},data.blocks.lenDist(1:end-1));
       blockLen_all = cat(2,blockLen_all,data.blocks.lenDist(1:end-1));
       end
        catch
            keyboard
        end
    end
    try
    blockNum_avg(na) = mean(blockNum{na});
    blockNum_sem(na) = nanstd(blockNum{na})./sqrt(numel(blockNum{na})); 
    blockNum_std(na) = nanstd(blockNum{na});
    catch
        
    blockNum_avg(na) = NaN;
    blockNum_sem(na) = NaN;
    blockNum_std(na) = NaN;
    end
     
  %  blockLen_avg(na) = mean(blockLen{na}); 
    
end
end
