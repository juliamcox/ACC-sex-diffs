clear all

behavfilename = 'P1E628_LON.mat';
syncfilename = 'P1_06282017_2.mat'; 


filebasename = dir('*downsamp*');
filebasename = filebasename.name;
filebasename = filebasename(1:strfind(filebasename,'.')-1);


frameRate = input('frame rate');
nFramesLog = input('n frames log');
droppedCount = input('dropped count');
droppedIdx = input('dropped idx');
ledPower = input('led power');
gain = input('gain');
startRec = input('startRec');
stopRec = input('stopRec');
recDur = input('recDur');
exposureTime = input('exposureTime');
save('info.mat');
