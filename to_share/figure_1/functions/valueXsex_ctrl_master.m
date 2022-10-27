clear all
close all
%% Parameters
cohort          = 'ACC_DMS_nphr';
aids_f = generateAnimalList('ACC_DMS_nphr_female');
aids_f = cat(1,aids_f,generateAnimalList('ACC_DMS_nphr_yfp_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d1_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_nphr_d2_female'));
aids_f = cat(1,aids_f,generateAnimalList('DMS_yfp_female')); 
%aids_f = cat(1,aids_f,generateAnimalList('imaging_female'));

aids_m = generateAnimalList('ACC_DMS_nphr_male');
aids_m = cat(1,aids_m,generateAnimalList('ACC_DMS_nphr_yfp_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d1_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_nphr_d2_male'));
aids_m = cat(1,aids_m,generateAnimalList('DMS_yfp_male')); 
%aids_m = cat(1,aids_m,generateAnimalList('imaging_male'));

cohort = 'all';

ext             = 'TAB'; 
sessionLength   = 'long';

binNum          = 4;
binNum_choice   = 9;
zscoreFlag      = 0; 
valType         = {'qChosenDiff'};
valType_choice  = {'qDiff'};
latencyType     = {'trialStart'}; 
perfThresh      = 0.1; 
basefilename    = fullfile(whereAreWe('bucket'), 'Manuscript_figures','fig1');
cutoff          = inf; 
lowCutoff       = 0; 
qFile           = 'qLearn_session_all.mat';
bandwidth       = .5; % bandwidth for kernel density function
ptiles          = [20:20:100];%
fext            = 'fig1';
%% Which functions to run

extractFlag(1) = 0; % [(1)extract latencies and choice in value bins and by previous outcome (Figure 1E,F Supp figure 1A,F; 
plotFlag(1)    = 1; 

extractFlag(2) = 0; %(2)plot Q-model parameters (Figure 1D);
plotFlag(2)    =0;  

extractFlag(3) = 0; %3)plot latency distribution by outcome(supp fig 1B,C);
plotFlag(3)    = 0;   

extractFlag(4) = 0; %(4)fit inverse gaussian (supp fig 1D-E) ]
plotFlag(4) = 1;

extractFlag(5) = 0; % (5) 2D latency x chosen and unchosen value 
plotFlag(5)    = 0;


%% Extract data

if extractFlag(1)
     ctrlXvalueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext);
end

if extractFlag(2)
    plotQValueParameters(aids_f,aids_m,sessionLength,perfThresh,qFile,ext);
    plotFlag(2) = 0; %function plots already
end

if extractFlag(3) || extractFlag(4)
    allCtrlData(aids_f,aids_m,sessionLength,perfThresh,qFile,basefilename)
end

if extractFlag(5)
    chosenUnchosen_2d(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext)
end
%% Plot latency x value x sex 

if plotFlag(1)
    % try to load extracted value, if not saved, extract
    try
        load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    catch
        ctrlXvalueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,'fig1');        
        load(fullfile(basefilename,['ctrlLatency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
        load(fullfile(basefilename,['ctrlLatency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_perfThresh_' num2str(perfThresh) '_cutoff_' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']))
    end
    stats_latency = latencyXvalueXsex_plot(latency_f, latency_m, valType,latencyType,zscoreFlag,bins);
    stats_latency_long = latencyXvalueXsex_plot_longTrials(latency_f, latency_m, valType,latencyType,zscoreFlag,bins);
    
    try
        load(fullfile(basefilename,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
        load(fullfile(basefilename,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
    catch
        ctrlXvalueXsex(zscoreFlag,binNum,binNum_choice, cohort,sessionLength,perfThresh,cutoff,lowCutoff,qFile,aids_m,aids_f,fext);
        load(fullfile(basefilename,['ctrlChoice_f_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
        load(fullfile(basefilename,['ctrlChoice_m_' cohort '_binNum' num2str(binNum_choice) '_perfThresh_' num2str(perfThresh) '.mat']))
    end
    stats_choice = choiceXvalueXsex_plot(choice_f, choice_m, valType_choice,binsChoice);
end

if plotFlag(2)
    plotQValueParameters(aids_f,aids_m,sessionLength,perfThresh,qFile,ext);
end

if plotFlag(3) 
    stats_hist=logLatencyHistograms_ctrl(qFile,sessionLength,perfThresh,aids_f,aids_m,basefilename);
end

if plotFlag(4)
    stats_invGauss= fitLatencyDistributions_inverseGaussian_control(aids_f,aids_m,sessionLength,perfThresh,qFile,basefilename,'none');
end
