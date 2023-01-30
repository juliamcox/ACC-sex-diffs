function fbasename = whereAreWe(loc)

loc = lower(loc);

switch loc
    case 'basis_sets'
        if ismac
            fbasename = [];
            error('Add directory location for data (in general_code folder)')
        elseif isunix
            fbasename = [];
            error('Add directory location for data (in general_code folder)')
        elseif ispc
            fbasename = [];
            error('Add directory location for data (in general_code folder)')
        end
    case 'data'
        if ismac
            fbasename = [];
            error('Add directory location for data')
        elseif isunix
            fbasename = [];
            error('Add directory location for data')
        elseif ispc
            fbasename = [];
            error('Add directory location for data')
        end
    case 'behavior'
        if ismac
            fbasename = [];
            error('Add directory location for data (Operant folder)')
        elseif isunix
            fbasename = [];
            error('Add directory location for data (Operant folder)')
        elseif ispc
            fbasename = [];
            error('Add directory location for data (Operant folder)')
        end
    case 'imaging'
        if ismac
            fbasename = [];
            error('Add directory location for data (ACC_DMS_imaging folder)')
        elseif isunix
            fbasename = [];
            error('Add directory location for data (ACC_DMS_imaging folder)')
        elseif ispc
            fbasename = [];
            error('Add directory location for data (ACC_DMS_imaging folder)')
        end
    case 'bilinear'
        if ismac
            fbasename = [];
            error('Add directory location for data (in ACC_DMS_imaging\bilinearRegression)')
        elseif isunix
            fbasename = [];
            error('Add directory location for data (in ACC_DMS_imaging\bilinearRegression)')
        elseif ispc
            fbasename = [];
            error('Add directory location for data (in ACC_DMS_imaging\bilinearRegression)')
        end
        
  
end
