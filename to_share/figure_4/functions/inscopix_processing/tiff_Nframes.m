function nZ = tiff_Nframes(filename)
%
% n = tiff_frames(filename)
%
% Returns the number of slices in a TIFF stack.
%
% copied from cellSort toolbox

status = 1; nZ=0;
jstep = 10^3;
while status
    try
        nZ=nZ+jstep;
        imread(filename,nZ);
    catch
        if jstep>1
            nZ=nZ-jstep;
            jstep = jstep/10;
        else
            nZ=nZ-1;
            status = 0;
        end
    end
end
end