function fbasename = whereAreWe(loc)

loc = lower(loc); 

switch loc
    case 'basis_sets'
        if ismac
            fbasename = '/Users/julia/Documents/code/ACC-sex-diffs/functions/general_code/basis_sets/';
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
    case 'data'
        if ismac
            fbasename = '/Users/julia/Documents/code/ACC-sex-diffs/processed_data';
        elseif isunix
            fbasename = [];
            error('Add directory location for data') 
        elseif ispc
            fbasename = [];
            error('Add directory location for data') 
        end
    case 'behavior'
        if ismac
            % fbasename = '/Volumes/witten/Julia/Operant/';
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/processed_data/Operant';
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
    case 'imaging'
        if ismac
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/processed_data/ACC_DMS_imaging/';
            %  fbasename = '/Volumes/witten/Julia/DMS_Bandit';
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
    case 'bilinear'
        if ismac
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/processed_data/ACC_DMS_imaging/bilinearRegression/';
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
        
        %     case 'human_bandit'
%         if ismac
%             fbasename = '/Volumes/witten/Julia/human_bandit/';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Julia/human_bandit/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Julia\human_bandit\';
%         end
%     case 'figurecode'
%         if ismac
%             fbasename = '/Users/julia/Documents/code/ACC-sex-diffs/';
%            % fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Paper_Code/Cox_2022/to_share/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Paper_Code\Cox_2022\to_share\';
%         end
%     case 'imaging'
%         if ismac
%           fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/processed_data/ACC_DMS_imaging/';
%        %  fbasename = '/Volumes/witten/Julia/DMS_Bandit';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Paper_Code/Cox_2022/to_share/processed_data/ACC_DMS_imaging/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Paper_Code\Cox_2022\to_share\processed_data\ACC_DMS_imaging/';
%         end
%     case 'behavior'
%         if ismac
%            % fbasename = '/Volumes/witten/Julia/Operant/';
%            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/processed_data/Operant';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Julia/Operant/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Julia\Operant\';
%         end
%     case 'bilinear'
%         if ismac
%             fbasename = '/Volumes/witten/Julia/DMS_Bandit/ACC_DMS_imaging/bilinearRegression/';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Julia/DMS_Bandit/ACC_DMS_imaging/bilinearRegression/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Julia\DMS_Bandit\ACC_DMS_imaging\bilinearRegression\';
%         end
%        
%     case 'bucket'
%         if ismac
%             fbasename = '/Volumes/witten/Julia/';
%         elseif isunix
%             fbasename = '/mnt/cup/labs/witten/Julia/';
%         elseif ispc
%             fbasename = '\\cup.pni.princeton.edu\witten\Julia\';
%         end
  
end
