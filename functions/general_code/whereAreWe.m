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
            fbasename ='/Volumes/witten/Paper_Code/Cox_2022/to_share/data';
        elseif isunix
            fbasename = [];
            error('Add directory location for data') 
        elseif ispc
            fbasename = [];
            error('Add directory location for data') 
        end
%     case 'behavior'
%         if ismac
%             % fbasename = '/Volumes/witten/Julia/Operant/';
%             fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/data/Operant';
%         elseif isunix
%             fbasename = [];
%             error('Add directory location for data')
%         elseif ispc
%             fbasename = [];
%             error('Add directory location for data')
%         end
    case 'imaging'
        if ismac
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/data/ACC_DMS_imaging/';
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
            fbasename = '/Volumes/witten/Paper_Code/Cox_2022/to_share/data/ACC_DMS_imaging/bilinearRegression/';
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
     
  
end
