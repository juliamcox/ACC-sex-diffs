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
    case 'imaging'
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
    case 'bilinear'
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
        
  
end
