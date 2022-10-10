function cnmfeDataConvert(dataLoc,recs)



%% extract and save metadata from .xml file 


% downsampled file 

filename = dir(fullfile(dataLoc,'*downsamp.tif*'));

idx = strfind(filename.name, '.');
filename = filename.name(1:idx-1);
fprintf([filename '\n']);
fprintf('\ncompiling metadata...\n')


if ~exist(fullfile(dataLoc,'info.mat'))
try
    % retrieve some info from inscopix-generated log file
    fid = fopen(fullfile(dataLoc, [filename(1:end-9) '.xml']));
    temp = textscan(fid,'%q');
catch
    fid = fopen(fullfile(dataLoc, 'session.json'));
    temp = textscan(fid,'%q');
end

% frame rate
colidx = (strfind(temp{1},'fps'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
frameRate = str2num(tempX(divs(1)+1:divs(end)-1));
% number of frames
colidx = (strfind(temp{1},'frames'));
rowidx = find(cellfun('isempty', colidx)==0);
rowidx = rowidx(1);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
nframesLog = str2num(tempX(divs(1)+1:divs(end)-1));
% dropped frames
colidx = (strfind(temp{1},'dropped_count'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
droppedCount = str2num(tempX(divs(1)+1:divs(end)-1));

if droppedCount ~= 0    
    colidx = (strfind(temp{1},'dropped'));
    rowidx = find(cellfun('isempty', colidx)==0);
    rowidx = rowidx(2);
    rowidxEnd = find(cellfun('isempty', strfind(temp{1}(rowidx:end), '<')) == 0)+rowidx-1;
    tempX = temp{1}(rowidx:rowidxEnd);
    tempX = cell2mat(tempX');
    divs = [strfind(tempX, '>') strfind(tempX, '<')];
    
    droppedIdx = str2num(tempX(divs(1)+1:divs(end)-1));
    if length(droppedIdx) ~= droppedCount
        keyboard
    end
elseif droppedCount == 0
    droppedIdx = []; 
end
% exposure
colidx = (strfind(temp{1},'exposure'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
exposureTime = str2num(tempX(divs(1)+1:divs(end)-1));
% LED power
colidx = (strfind(temp{1},'led_power'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
ledPower = str2num(tempX(divs(1)+1:divs(end)-1));
% gain
colidx = (strfind(temp{1},'gain'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = temp{1}{rowidx};
divs = [strfind(tempX, '>') strfind(tempX, '<')];
gain = str2num(tempX(divs(1)+1:divs(end)-1));
% start recording
colidx = (strfind(temp{1},'record_start'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = [temp{1}{rowidx+3} ' ' temp{1}{rowidx+4}(1:2)];
startRec = vpa(datenum(tempX));
% stop recording
colidx = (strfind(temp{1},'record_end'));
rowidx = find(cellfun('isempty', colidx)==0);
tempX = [temp{1}{rowidx+3} ' ' temp{1}{rowidx+4}(1:2)];
stopRec = vpa(datenum(tempX));
recDur = (stopRec-startRec)*24*3600; %in sec

fclose(fid);
clear temp rowidx colidx tempX fid
save(fullfile(dataLoc, 'info.mat')) % save metadata
end

%% check to see if the data has previously been saved in imstack.mat

if exist(fullfile(dataLoc, 'imstack.mat')) == 0
    % convert tif stack
    fprintf('Converting tif stack to mat...\n')
    
    xx = strfind(recs, '\');
    if isempty(xx)
        xx = strfind(recs,'/');
    end
    imstack = loadTiffStack(dataLoc, []);
    Y = double(imstack); clear imstack
    
    
else
    fprintf('Converting imstack...\n')
    imstackObj = matfile(fullfile(dataLoc,'imstack.mat'));
    imstackObj.Properties.Writable = true;
    if length(whos(imstackObj)) > 1
        
        filelist = dir([dataLoc recs '/*downsamp*.tif']);
        count = 1;
        Y = [];
        x = length(whos(imstackObj));
        for i = 1:length(whos(imstackObj))
            fprintf('%d/%d....\n', i, length(whos(imstackObj)));
            eval(sprintf('Y = cat(3, Y, double(imstackObj.imstack_part%d));', i));
        end
    else
        Y = double(imstackObj.imstack);
    end
    
    
    
end

info = matfile(fullfile(dataLoc,'info.mat'));
Fs = info.frameRate;             % frame rate
Ysiz = size(Y);
save(fullfile(dataLoc,'Y.mat'), 'Y', 'Ysiz', 'Fs', '-v7.3')


