function extractBlockLength(aids_f,aids_m,perfThresh,qFile,sessionLength)

[blockNum,blockNum_avg,blockNum_sem,blockNum_std,blockLen]=extractSession(aids_f,perfThresh,qFile,sessionLength);
[blockNum_m,blockNum_avg_m,blockNum_sem_m,blockNum_std_m,blockLen_m]=extractSession(aids_m,perfThresh,qFile,sessionLength);

keyboard
end


function [blockNum,blockNum_avg,blockNum_sem,blockNum_std,blockLen_all,blockLen]=extractSession(aids,perfThresh,qFile,sessionLength)
basefilename = fullfile(whereAreWe('bucket'), 'Operant');
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
