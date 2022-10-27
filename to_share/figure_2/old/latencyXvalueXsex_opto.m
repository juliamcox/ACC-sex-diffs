function latencyXvalueXsex_opto(epochs,ext,zscoreFlag,cohort,binNum,cutoff,lowCutoff,qFile)

%%% Figure 4 %%%

% Extract trial initiation time, lever press latency or nose poke withdrawal as a function of Q-value

%%% Inputs
% ext: file extension (LAS,UNI et)
% epochs: laser conditions to plot
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% zscoreFlag: 0 or 1 to zscore resposnse times
% cohort: basename for cohort (e.g. 'ACC_DMS_nphr') for animal lists generated with generateAnimalList

%%% Dependencies
% Stan Q model converted to matlab with qValue_individual_fromMat.mat
% generateAnimalList: to create subject lists
% whereAreWe: to create file path
% latency_value_opto: embedded function to extract latency information 
% getIndices
% responseTimesGUI
% concatValue extract sessions 

%% Parameters 

savehere = fullfile(whereAreWe('bucket'), 'Manuscript_figures','fig4');

% Generate male, female, opsin, yfp subject  lists 
aids_opsin_m = generateAnimalList([cohort '_male']);
aids_yfp_m = generateAnimalList([cohort '_yfp_male']);

aids_opsin_f = generateAnimalList([cohort '_female']);
aids_yfp_f = generateAnimalList([cohort '_yfp_female']);
% Load combined D1/D2 YFP group if DMS 
if ~contains(cohort,'ACC')
    aids_yfp_f = generateAnimalList('DMS_yfp_female');
    aids_yfp_m = generateAnimalList('DMS_yfp_male');
    
end

ext_opsin_f = repmat({ext},size(aids_opsin_f));
ext_yfp_f = repmat({ext},size(aids_yfp_f));

ext_opsin_m = repmat({ext},size(aids_opsin_m));
ext_yfp_m = repmat({ext},size(aids_yfp_m));
%% Extract and save latency x value 

[latency_m,bins] = latency_value_opto(aids_opsin_m,epochs,ext_opsin_m,zscoreFlag,qFile,binNum,cutoff,lowCutoff);
save(fullfile(savehere,['latency_m_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_cutoff' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_m','bins')
[latency_m_yfp,bins] = latency_value_opto(aids_yfp_m,epochs,ext_yfp_m,zscoreFlag,qFile,binNum,cutoff,lowCutoff);
save(fullfile(savehere,['latency_m_yfp_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_cutoff' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_m_yfp','bins')
[latency_f,bins] = latency_value_opto(aids_opsin_f,epochs,ext_opsin_f,zscoreFlag,qFile,binNum,cutoff,lowCutoff);
save(fullfile(savehere,['latency_f_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_cutoff' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_f','bins')
[latency_f_yfp,bins] = latency_value_opto(aids_yfp_f,epochs,ext_yfp_f,zscoreFlag,qFile,binNum,cutoff,lowCutoff);
save(fullfile(savehere,['latency_f_yfp_zscore' num2str(zscoreFlag) '_' cohort '_binNum' num2str(binNum) '_cutoff' num2str(cutoff) '_lowCutoff' num2str(lowCutoff) '.mat']), 'latency_f_yfp','bins')

end
%% Extract trial initiation latency X value 
function [latency,bins] = latency_value_opto(aids,epochs,ext,zscoreFlag,qFile,binNum,cutoff,lowCutoff)

%aids: list of animal IDs
%ext: file extension (LAS,UNI et) 
%epochs: laser conditions to extract 
%zscoreFlag: 0 or 1 to zscore resposnse times 
%qFile: which q values to extract
%binNum: binNum for value
%cutoff: exclude long ITIs?

bins = [];


basefilename = fullfile(whereAreWe('bucket'), 'Operant');

% Value types
valTypes = {'qChosenDiff_norm';'qChosenDiff_abs';'qDiff';'qTot';'qChosen';'qUnchosen';'qChosenDiff';'qIpsiDiff';'qDiff_prev';'qTot_prev';'qChosen_prev';'qUnchosen_prev';'qChosenDiff_prev';'qIpsiDiff_prev'}; 
latencyTypes = {'trialStart';'withdraw';'leverPress'}; 

% Initialize variables for response times X value
for nv = 1:numel(valTypes)
    for nl = 1:numel(latencyTypes)
        eval(sprintf('latency.%s_%s = nan(numel(aids),binNum);',latencyTypes{nl}, valTypes{nv}));
        eval(sprintf('latency.%s_%s_count = nan(numel(aids),binNum);',latencyTypes{nl}, valTypes{nv}));
        for ne = 1:numel(epochs)
            eval(sprintf('latency.%s_%s_%s = nan(numel(aids),binNum);',latencyTypes{nl}, valTypes{nv}, epochs{ne}));
            eval(sprintf('latency.%s_%s_%s_count = nan(numel(aids),binNum);',latencyTypes{nl}, valTypes{nv}, epochs{ne}));
        end
    end
end

for na = 1:numel(aids)
    
    load(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat'])); % Load behavior data
    
    % Find laser trials
    if ~isfield(data, 'idx')
        data = getIndices(data);
        save(fullfile(basefilename, aids{na}, ['/laserSummary_' ext{na} '.mat']), 'data')
    end
    
    % If response times have not been calculated 
    if ~isfield(data, 'rt')
        data = responseTimesGUI(data);
    end
    if zscoreFlag == 1
        trialStart = data.rt.nosePoke_zscore;
        leverPress = data.rt.leverPress_zscore;
        withdraw   = data.rt.withdraw_zscore;
        
        withdraw(data.rt.withdraw > cutoff | data.rt.withdraw < lowCutoff) = NaN; %cutoff in seconds
        trialStart(data.rt.nosePoke > cutoff | data.rt.nosePoke < lowCutoff) = NaN;
        leverPress(data.rt.leverPress > cutoff | data.rt.leverPress < lowCutoff) = NaN;
        
    elseif zscoreFlag == 2     
        data.rt.nosePoke(data.rt.nosePoke==0) = .01;
        trialStart = log(data.rt.nosePoke);
        leverPress = log(data.rt.leverPress);
        withdraw   = log(data.rt.withdraw);
        
        trialStart(data.rt.nosePoke > cutoff | data.rt.nosePoke < lowCutoff) = NaN;
        leverPress(data.rt.leverPress > cutoff | data.rt.leverPress < lowCutoff) = NaN;
        withdraw(data.rt.withdraw > cutoff | data.rt.withdraw < lowCutoff) = NaN;
    else
        trialStart = data.rt.nosePoke;
        leverPress = data.rt.leverPress;
        withdraw   = data.rt.withdraw;
        trialStart(data.rt.nosePoke > cutoff | data.rt.nosePoke < lowCutoff) = NaN;
        leverPress(data.rt.leverPress > cutoff | data.rt.leverPress < lowCutoff) = NaN;
        withdraw(data.rt.withdraw > cutoff | data.rt.withdraw < lowCutoff) = NaN;
     end
%     
% % %     
%     if ~isfield(data,'QIpsi')
%         % Load Q-values and extract laser sessions
%         load(fullfile(basefilename, aids{na}, qFile));
%         data = concatValue(data,qLearn,aids{na});
%         data.qFile = qFile;
%         save(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat']), 'data')
%     elseif ~contains(data.qFile,qFile)
%        Load Q-values and extract laser sessions
        load(fullfile(basefilename, aids{na}, qFile));
        data = concatValue(data,qLearn,aids{na});
        data.qFile = qFile;
        save(fullfile(basefilename, aids{na}, ['laserSummary_' ext{na} '.mat']), 'data')
 %  end
    
    qDiff = data.QRight-data.QLeft; % relative value
    qTot =  data.QRight+data.QLeft; % total value
    clear qChosen qUnchosen qChosenDiff
    qChosen = nan(size(data.choice));
    qUnchosen = nan(size(data.choice));
    qChosen(data.choice==1) = data.QLeft(data.choice==1);
    qChosen(data.choice==0) = data.QRight(data.choice==0);
    qUnchosen(data.choice==1) = data.QRight(data.choice==1);
    qUnchosen(data.choice==0) = data.QLeft(data.choice==0);
    qChosenDiff_abs = abs(qChosen-qUnchosen);
    qChosenDiff_norm = (qChosen-qUnchosen)./(qChosen+qUnchosen);
    qChosenDiff = (qChosen-qUnchosen);
    qIpsiDiff = data.QIpsi-data.QContra;
    qChosen_prev = [NaN qChosen(1:end-1)];
    qTot_prev = [NaN qDiff(1:end-1)];
    qChosenDiff_prev = [NaN qChosenDiff(1:end-1)];
    qDiff_prev = [NaN qDiff(1:end-1)];
    qUnchosen_prev = [NaN qUnchosen(1:end-1)];
    qIpsiDiff_prev = [NaN qIpsiDiff(1:end-1)];

    [~,bins.qDiff,qDiff_bins] = histcounts(qDiff,linspace(-1,1,binNum+1));
    [~,bins.qTot,qTot_bins] = histcounts(qTot,linspace(0,2,binNum+1));
    [~,bins.qChosen,qChosen_bins] = histcounts(qChosen,linspace(0,1,binNum+1));
    [~,bins.qUnchosen,qUnchosen_bins] = histcounts(qUnchosen,linspace(0,1,binNum+1));
    [~,bins.qChosenDiff,qChosenDiff_bins] = histcounts(qChosenDiff,linspace(-1,1,binNum+1));
    [~,bins.qIpsiDiff,qIpsiDiff_bins] = histcounts(qIpsiDiff,linspace(-1,1,binNum+1));
    [~,bins.qChosenDiff_abs,qChosenDiff_abs_bins] = histcounts(qChosenDiff_abs,linspace(0,1,binNum+1));
    [~,bins.qChosenDiff_norm,qChosenDiff_norm_bins] = histcounts(qChosenDiff_norm,linspace(-1,1,binNum+1));

    
    [~,bins.qDiff_prev,qDiff_prev_bins] = histcounts(qDiff_prev,linspace(-1,1,binNum+1));
    [~,bins.qTot_prev,qTot_prev_bins] = histcounts(qTot_prev,linspace(0,2,binNum+1));
    [~,bins.qChosen_prev,qChosen_prev_bins] = histcounts(qChosen_prev,linspace(0,1,binNum+1));
    [~,bins.qUnchosen_prev,qUnchosen_prev_bins] = histcounts(qUnchosen_prev,linspace(0,1,binNum+1));
    [~,bins.qChosenDiff_prev,qChosenDiff_prev_bins] = histcounts(qChosenDiff_prev,linspace(-1,1,binNum+1));
    [~,bins.qIpsiDiff_prev,qIpsiDiff_prev_bins] = histcounts(qIpsiDiff_prev,linspace(-1,1,binNum+1));
    
    % Bin values in quantiles by session 
%     clear qChosen_quant qUnchosen_quant qDiff_quant qChosenDiff_quant qTot_quant
%     breaks = find(data.choice==-10);
%     session_start = [1 breaks(2:2:end-1)+1];
%     session_end = [breaks(1:2:end-1)-1 numel(data.choice)];
%     for ns = 1:numel(session_start)
%         if sum(isnan(qChosen(session_start(ns):session_end(ns)))) == numel(session_start(ns):session_end(ns))
%             qChosen_quant{ns} = nan(size(session_start(ns):session_end(ns)));
%             qUnchosen_quant{ns} = nan(size(session_start(ns):session_end(ns)));
%             qDiff_quant{ns} = nan(size(session_start(ns):session_end(ns)));
%             qChosenDiff_quant{ns} = nan(size(session_start(ns):session_end(ns)));
%             qTot_quant{ns} = nan(size(session_start(ns):session_end(ns)));
%         else
%             tempBins = prctile(qChosen(session_start(ns):session_end(ns)), [0 100/binNum:100/binNum:100]);
%             [~,~,qChosen_quant{ns}] = histcounts(qChosen(session_start(ns):session_end(ns)),tempBins);
%             
%             tempBins = prctile(qUnchosen(session_start(ns):session_end(ns)), [0 100/binNum:100/binNum:100]);
%             [~,~,qUnchosen_quant{ns}] = histcounts(qUnchosen(session_start(ns):session_end(ns)),tempBins);
%             tempBins = prctile(qTot(session_start(ns):session_end(ns)), [0 100/binNum:100/binNum:100]);
% try
%             [~,~,qTot_quant{ns}] = histcounts(qTot(session_start(ns):session_end(ns)),tempBins);
% catch
%     keyboard
% end
%             tempBins = prctile(qChosenDiff(session_start(ns):session_end(ns)), [0 100/binNum:100/binNum:100]);
%             [~,~,qChosenDiff_quant{ns}] = histcounts(qChosenDiff(session_start(ns):session_end(ns)),tempBins);
%             tempBins = prctile(qDiff(session_start(ns):session_end(ns)), [0 100/binNum:100/binNum:100]);
%             [~,~,qDiff_quant{ns}] = histcounts(qDiff(session_start(ns):session_end(ns)),tempBins);
%             
%         end
%     end
%     qChosen_quant = cell2mat(qChosen_quant);
%     qUnchosen_quant = cell2mat(qUnchosen_quant);
%     qDiff_quant = cell2mat(qDiff_quant);
%     qChosenDiff_quant = cell2mat(qChosenDiff_quant);
%     qTot_quant = cell2mat(qTot_quant);

    
%       % Bin values in quantiles
%     
    tempBins = [prctile(qChosen,[0 100/binNum:100/binNum:100])];
    [~,~,qChosen_quant] = histcounts(qChosen,tempBins);
    
    tempBins = [prctile(qUnchosen,[0 100/binNum:100/binNum:100])];
    [~,~,qUnchosen_quant] = histcounts(qUnchosen,tempBins);
    try
    tempBins = [prctile(qTot,[0 100/binNum:100/binNum:100])];
    [~,~,qTot_quant] = histcounts(qTot,tempBins);
    catch
        keyboard
    end
    
    tempBins = [prctile(qDiff,[0 100/binNum:100/binNum:100])];
    [~,~,qDiff_quant] = histcounts(qDiff,tempBins);
    
    tempBins = [prctile(qChosenDiff,[0 100/binNum:100/binNum:100])];
    [~,~,qChosenDiff_quant] = histcounts(qChosenDiff,tempBins);
     tempBins = [prctile(qChosenDiff_abs,[0 100/binNum:100/binNum:100])];
    [~,~,qChosenDiff_abs_quant] = histcounts(qChosenDiff_abs,tempBins);
    
     tempBins = [prctile(qChosenDiff_norm,[0 100/binNum:100/binNum:100])];
    [~,~,qChosenDiff_norm_quant] = histcounts(qChosenDiff_norm,tempBins);
    
    tempBins = [prctile(qIpsiDiff,[0 100/binNum:100/binNum:100])];
    [~,~,qIpsiDiff_quant] = histcounts(qIpsiDiff,tempBins);
%     


     tempBins = [prctile(qChosen_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qChosen_prev_quant] = histcounts(qChosen_prev,tempBins);
    
    tempBins = [prctile(qUnchosen_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qUnchosen_prev_quant] = histcounts(qUnchosen_prev,tempBins);
    try
    tempBins = [prctile(qTot_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qTot_prev_quant] = histcounts(qTot_prev,tempBins);
    catch
        keyboard
    end
    
    tempBins = [prctile(qDiff_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qDiff_prev_quant] = histcounts(qDiff_prev,tempBins);
    
    tempBins = [prctile(qChosenDiff_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qChosenDiff_prev_quant] = histcounts(qChosenDiff_prev,tempBins);
    
    
    tempBins = [prctile(qIpsiDiff_prev,[0 100/binNum:100/binNum:100])];
    [~,~,qIpsiDiff_prev_quant] = histcounts(qIpsiDiff_prev,tempBins);
    
    % Divide response times by current trial value for laser and non-laser trials
    for nv = 1:numel(valTypes)
        thisbins = eval(sprintf('bins.%s;',valTypes{nv}));
        for nb = 1:numel(thisbins)
            for nl = 1:numel(latencyTypes)
                % non-laser trials
                idx = data.trialIdx;
                idx = idx(ismember(idx-1,data.trialIdx));
                eval(sprintf('idx = idx(ismember(idx,find(%s_bins == nb)));',valTypes{nv}));
                eval(sprintf('latency.%s_%s(na,nb) = nanmean(%s(idx));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                %eval(sprintf('latency.%s_%s(na,nb) = prctile(%s(idx),80);', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                eval(sprintf('latency.%s_%s_count(na,nb) = numel(idx);', latencyTypes{nl}, valTypes{nv}));
                % non-laser trials
                idx = data.trialIdx;
                idx = idx(ismember(idx-1,data.trialIdx));
                eval(sprintf('idx = idx(ismember(idx,find(%s_quant == nb)));',valTypes{nv}));
                eval(sprintf('latency.%s_%s_quant(na,nb) = nanmean(%s(idx));', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                %eval(sprintf('latency.%s_%s_quant(na,nb) = prctile(%s(idx),80);', latencyTypes{nl}, valTypes{nv}, latencyTypes{nl}));
                eval(sprintf('latency.%s_%s_quantCount(na,nb) = numel(idx);', latencyTypes{nl}, valTypes{nv}));


                
                % laser trials
                for ne = 1:numel(epochs)
                    idx_laser = eval(sprintf('data.trialIdx_laser%s',epochs{ne}));
                    eval(sprintf('idx_laser = idx_laser(ismember(idx_laser,find(%s_bins == nb)));',valTypes{nv}));
                    eval(sprintf('latency.%s_%s_%s(na,nb) = nanmean(%s(idx_laser));',latencyTypes{nl}, valTypes{nv}, epochs{ne},latencyTypes{nl}));
                    %eval(sprintf('latency.%s_%s_%s(na,nb) = prctile(%s(idx_laser),80);',latencyTypes{nl}, valTypes{nv}, epochs{ne},latencyTypes{nl}));
                    eval(sprintf('latency.%s_%s_%s_count(na,nb) = numel(idx_laser);',latencyTypes{nl}, valTypes{nv}, epochs{ne}));
                    
                    idx_laser = eval(sprintf('data.trialIdx_laser%s',epochs{ne}));
                    eval(sprintf('idx_laser = idx_laser(ismember(idx_laser,find(%s_quant == nb)));',valTypes{nv}));
                    eval(sprintf('latency.%s_%s_%s_quant(na,nb) = nanmean(%s(idx_laser));', latencyTypes{nl}, valTypes{nv}, epochs{ne}, latencyTypes{nl}));
                    %eval(sprintf('latency.%s_%s_%s_quant(na,nb) = prctile(%s(idx_laser),80);', latencyTypes{nl}, valTypes{nv}, epochs{ne}, latencyTypes{nl}));
                    eval(sprintf('latency.%s_%s_%s_quantCount(na,nb) = numel(idx_laser);', latencyTypes{nl}, valTypes{nv}, epochs{ne}));

                end
            end
        end
    end
    
    
    %% Long trials (if needed, needs to be updated for change in cutoff indexing 
%     for nl = 1:numel(latencyTypes)
%         % non-laser trials
%         eval(sprintf('idx_long = data.trialIdx(%s(data.trialIdx)>cutoff);',latencyTypes{nl}));
%         eval(sprintf('latency.%s_longPer(na) = numel(idx_long)/numel(data.trialIdx);',latencyTypes{nl}));
%         for ne = 1:numel(epochs)
%             idx_laser = eval(sprintf('data.trialIdx_laser%s',epochs{ne}));
%             eval(sprintf('idx_laser_long = idx_laser(%s(idx_laser)>cutoff);',latencyTypes{nl}))
%             eval(sprintf('latency.%s_%s_longPer(na) = numel(idx_laser_long)/numel(idx_laser);',latencyTypes{nl},epochs{ne})); 
%         end
%         
%     end

end
end