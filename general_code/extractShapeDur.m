function numShapeSess = extractShapeDur(aids)

basefilename = fullfile(whereAreWe('behavior'));

for na = 1:numel(aids)
    cd(fullfile(basefilename,aids{na}));
    fnames = dir('*.SHA');
    numShapeSess(na) = numel(fnames)-1;
end
    


end
