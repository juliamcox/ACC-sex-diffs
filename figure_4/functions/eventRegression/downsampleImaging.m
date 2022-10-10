function data_ds = downsampleImaging(frameRate_orig, frameRate_new,data)

if size(data,2) > size(data,1)
    data = data';
end

% downsample frame times (will not do anything if already in desired sampling frequency)
endPoint = length(data)/frameRate_orig;
currTaxis = linspace(0,endPoint,length(data))';
newTaxis = (0:1/frameRate_new:endPoint)';
if newTaxis(end)<endPoint; newTaxis = [newTaxis; endPoint]; end
data_ds = interp1(currTaxis,data,newTaxis); 
