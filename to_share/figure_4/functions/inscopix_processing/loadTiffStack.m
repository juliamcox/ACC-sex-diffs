function imstack = loadTiffStack(filebasename,frameRange,singFileFlag)

% imstack = loadTiffStack(filebasename,frameRange)
% loads tif stack from 1 or multiple files
% LP may 2013

if nargin == 1
    frameRange = [];
    singFileFlag = 0;
elseif nargin == 2
    singFileFlag = 0;
end
xx = strfind(filebasename, '/');
    if isempty(xx)
        strfind(filebasename, '\')
    end

% figure out how many frames per file and image size
if singFileFlag
    filelist{1} = filebasename;
else
    xx = strfind(filebasename, '/');
    if isempty(xx)
        strfind(filebasename, '\')
    end
    filelist = generateFilelist(fullfile(filebasename, [filebasename(xx(end)+1:end) '_downsamp']));
end


nZ = zeros(1,length(filelist));
%filelist = cellfun(@(x) x{1}, filelist, 'UniformOutput', false);
tic;
fprintf('loading tiff stack...\n')
for i = 1:length(filelist)
    nZ(i) = tiff_Nframes(fullfile(filebasename,filelist{i}));
end

temp = imread(fullfile(filebasename,filelist{i}),1);
[nX nY] = size(temp); clear temp

if isempty(frameRange); frameRange = [1 sum(nZ)]; end

if length(filelist) == 1
    %     imstack = TIFFStack(fullfile(filebasename, filelist{1}));
    %     imstack = imstack(:,:,frameRange(1):frameRange(end));
    
    FileTif= fullfile(filebasename, filelist{1});
    imstack=zeros(nX,nY,nZ,'uint16');
    
    TifLink = Tiff(FileTif, 'r');
    for i=1:nZ
        TifLink.setDirectory(i);
        imstack(:,:,i)=TifLink.read();
    end
    TifLink.close();
else
    
    imstack = zeros(nX,nY,frameRange(end));
    count = 1;
    for i = 1:length(filelist)
        
        FileTif= fullfile(filebasename, filelist{i});
        temp=zeros(nX,nY,nZ(i),'uint16');
        
        TifLink = Tiff(FileTif, 'r');
        for ii=1:nZ(i)
            TifLink.setDirectory(ii);
            temp(:,:,ii)=TifLink.read();
        end
        TifLink.close();
        
       % temp=TIFFStack(fullfile(filebasename,filelist{i}));
        if count+nZ(i)-1 <= frameRange(end)
            imstack(:,:,count:count+nZ(i)-1)=temp(:,:,:);
        else imstack(:,:,count:frameRange(end))=temp(:,:,1:frameRange(end)-count+1);
        end
        clear temp
        count = count + nZ(i);
    end
end
t = toc; tic
fprintf(' done after %1.2f minutes\n',t/60)
end


% generate downsampled file list
function filelist = generateFilelist(filebasename)
filelisttemp = dir([filebasename '*tif*']); count = 1;

for i = 1:size(filelisttemp,1)
    filelisttemp_temp{i,:} = filelisttemp(i).name;
    %     x = strfind(filelisttemp_temp(i,:),'f');
    %     filelisttemp_temp2(i,:) = filelisttemp_temp(i,1:x(end));
    
end
filelist = filelisttemp_temp;
% filelisttemp = filelisttemp_temp;
% for i = 1:size(filelisttemp,1)
%     if isempty(strfind(filelisttemp{i},'aligned'))
%         fl(count,:)=filelisttemp{i,:};
%         count = count + 1;
%     end
% end
% if size(fl,1) > 1
%     count = 1;
%     for i=[size(fl,1) 1:size(fl,1)-1] % invert because 1st file comes last on list
%         idx=strfind(fl{i,:},'f'); % last non-empty character of file name
%         filelist{count}=fl{i,1:idx(end)};
%         count = count +1;
%     end
% else
%     filelist{1}=fl{1,:};
% end
end
