
%% Parameters
sessionLength = 'long';
perfThresh    = 0.1;
binNum        = 9;
qFile         =  'qLearn_session_all_2022.mat';
aids_f        = {'I40';'N10';'D9';'T98';'D37';'D38';'D17';'I53';'I42';'I49';'I45';'N9';'T96';'I47';'T97';'I48'};
aids_m        = {'110';'N2';'105';'N5';'107';'D29';'I28';'109';'108';'N1';'D15';'D32';'T49';'I17';'104';'N4'};
%% Plot 
plotPsychCurves(sessionLength, perfThresh, qFile, binNum, aids_f, 1)
plotPsychCurves(sessionLength, perfThresh, qFile, binNum, aids_m, 0)