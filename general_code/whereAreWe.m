function fbasename = whereAreWe(loc)

loc = lower(loc); 

switch loc
    case 'human_bandit'
        if ismac
            fbasename = '/Volumes/witten/Julia/human_bandit/';
        elseif isunix
            fbasename = '/mnt/cup/labs/witten/Julia/human_bandit/';
        elseif ispc
            fbasename = '\\cup.pni.princeton.edu\witten\Julia\human_bandit\';
        end
    case 'figurecode'
        if ismac
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022';
        elseif isunix
            fbasename = '/mnt/cup/labs/witten/Paper_Code/Cox_2022';
        elseif ispc
            fbasename = '\\cup.pni.princeton.edu\witten\Paper_Code\Cox_2022\';
        end
    case 'imaging'
        if ismac
            fbasename = '/Volumes/witten/Julia/DMS_Bandit/ACC_DMS_imaging/';
        elseif isunix
            fbasename = '/mnt/cup/labs/witten/Julia/DMS_Bandit/ACC_DMS_imaging/';
        elseif ispc
            fbasename = '\\cup.pni.princeton.edu\witten\Julia\DMS_Bandit\ACC_DMS_imaging\';
        end
    case 'behavior'
        if ismac
            fbasename = '/Volumes/witten/Julia/Operant/';
        elseif isunix
            fbasename = '/mnt/cup/labs/witten/Julia/Operant/';
        elseif ispc
            fbasename = '\\cup.pni.princeton.edu\witten\Julia\Operant\';
        end
       
    case 'bucket'
        if ismac
            fbasename = '/Volumes/witten/Julia/';
        elseif isunix
            fbasename = '/mnt/cup/labs/witten/Julia/';
        elseif ispc
            fbasename = '\\cup.pni.princeton.edu\witten\Julia\';
        end
  
end
