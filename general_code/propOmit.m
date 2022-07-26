
% Parameters 
aids_f = generateAnimalList('ACC_DMS_nphr_female');
aids_f = cat(1,aids_f,generateAnimalList('ACC_DMS_nphr_yfp_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d1_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d2_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_yfp_female')); 
aids_f = cat(1,aids_f,generateAnimalList('imaging_female'));

aids_m = generateAnimalList('ACC_DMS_nphr_male');
aids_m = cat(1,aids_m,generateAnimalList('ACC_DMS_nphr_yfp_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d1_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d2_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_yfp_male')); 
aids_m = cat(1,aids_m,generateAnimalList('imaging_male')); 

ext             = 'TAB'; 
sessionLength  = 'both';

perfThresh      = 0.1; 
qFile           = 'qLearn_session_all_2022.mat';
basefilename = fullfile(whereAreWe('behavior'));

aids = cat(1,aids_f,aids_m);
%% Find proportion omitted trials 

for na = 1:numel(aids)
    try load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist; 
    catch
        valueExtraction_TAB(aids(na),qFile,sessionLength,perfThresh);
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist;
    end
    omitProp(na) = sum(omit)/sum(totalTrial);
end

% estrous mice 
qFile = 'qLearn_session_estrous_2022.mat';
sessionLength = 'long';
aids = generateAnimalList('estrous');
x = na+1;
for na = 1:numel(aids)
    try load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist; 
    catch
        valueExtraction_TAB(aids(na),qFile,sessionLength,perfThresh);
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        load(fullfile(basefilename,aids{na}, sprintf('valueTAB_flist_%s_perfThresh_%s_%s',sessionLength,num2str(perfThresh),qFile)))
        thisFlist{na} = flist;
    end
    omitProp(x) = sum(omit)/sum(totalTrial);
    x=x+1;
end