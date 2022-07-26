function [mmap, fmap] = maleFemaleColormap(maleC,femaleC,N)
if nargin < 3
    N = 100;
end


mmap = flipud([linspace(maleC(1),1,N); linspace(maleC(2),1,N); linspace(maleC(3),1,N)]');
fmap = flipud([linspace(femaleC(1),1,N); linspace(femaleC(2),1,N); linspace(femaleC(3),1,N)]');
