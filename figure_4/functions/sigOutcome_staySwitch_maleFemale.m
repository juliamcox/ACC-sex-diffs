function sigOutcome_staySwitch_maleFemale(rawFlag,frameRate,zscoreFlag,sigVer,sigEvent,qFile,events,eventDur,kernelVer)

%rawFlag: plot raw (1) or denoised (0) gcamp 
%zscoreFlag: zscore fluorescence? 
%sigVer: which version of regression to use (if only extracting significant neurons)
%sigEvent: which event 
%qFile: 
%events: what event to extract fluorescence (stay v. switch)
%eventDur: time before and after event to plot
%kernelVer: which regression version to plot kernel
%expt: which imaging cohort? 

if nargin < 10
    expt = 'ACC_DMS_imaging';
end

recs_f = generateAnimalList('ACC_DMS_imaging_female');
recs_m = generateAnimalList('ACC_DMS_imaging_male');

out_f = extractStaySwitch(recs_f,rawFlag,frameRate,zscoreFlag,sigVer,sigEvent,qFile,events,eventDur,kernelVer);
out_m = extractStaySwitch(recs_m,rawFlag,frameRate,zscoreFlag,sigVer,sigEvent,qFile,events,eventDur,kernelVer);
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

%plotExampleStaySwitch2(out_f,eventDur,plotParams.femaleC);

%% Plot   
% Plot mean activity for reward -> stay v switch and no reward -> stay v switch
for ne = 1:numel(events)
    % Plot mean activity across cells
    figure();
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out_f.dff_rew_switch{ne},2));
    subplot(1,2,1); hold on, box off
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willStay{ne}), std(out_f.dff_rew_willStay{ne})./sqrt(size(out_f.dff_rew_willStay{ne},1)),'LineProps',{'Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willSwitch{ne}), std(out_f.dff_rew_willSwitch{ne})./sqrt(size(out_f.dff_rew_willSwitch{ne},1)),'LineProps',{'--','Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willStay{ne}), std(out_m.dff_rew_willStay{ne})./sqrt(size(out_m.dff_rew_willStay{ne},1)),'LineProps',{'Color',plotParams.maleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willSwitch{ne}), std(out_m.dff_rew_willSwitch{ne})./sqrt(size(out_m.dff_rew_willSwitch{ne},1)),'LineProps',{'--','Color',plotParams.maleC})
    
    xlabel('Time from reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out_f.dff_rew_switch{ne},2));
    subplot(1,2,2); hold on, box off
    lh{1}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willStay{ne}), std(out_f.dff_nrew_willStay{ne})./sqrt(size(out_f.dff_nrew_willStay{ne},1)),'LineProps',{'Color',plotParams.femaleC});
    lh{2}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willSwitch{ne}), std(out_f.dff_nrew_willSwitch{ne})./sqrt(size(out_f.dff_nrew_willSwitch{ne},1)),'LineProps',{'--','Color',plotParams.femaleC});
     lh{3}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willStay{ne}), std(out_m.dff_nrew_willStay{ne})./sqrt(size(out_m.dff_nrew_willStay{ne},1)),'LineProps',{'Color',plotParams.maleC});
    lh{4}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willSwitch{ne}), std(out_m.dff_nrew_willSwitch{ne})./sqrt(size(out_m.dff_nrew_willSwitch{ne},1)),'LineProps',{'--','Color',plotParams.maleC});

    xlabel('Time from no reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    legend([lh{1}.mainLine lh{2}.mainLine lh{3}.mainLine lh{4}.mainLine],{'female: Stay', 'female:Switch','male: Stay', 'male:Switch'})
    
    
    %% Plot mean activity across cells (cs plus v csmius
    figure();
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out_f.dff_rew_switch{ne},2));
    subplot(2,2,1); hold on, box off
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willStay{ne}(logical(out_f.negKernel{ne}),:)), std(out_f.dff_rew_willStay{ne}(logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_rew_willStay{ne}(logical(out_f.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willSwitch{ne}(logical(out_f.negKernel{ne}),:)), std(out_f.dff_rew_willSwitch{ne}(logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_rew_willSwitch{ne}(logical(out_f.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willStay{ne}(logical(out_m.negKernel{ne}),:)), std(out_m.dff_rew_willStay{ne}(logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_rew_willStay{ne}(logical(out_m.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.maleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willSwitch{ne}(logical(out_m.negKernel{ne}),:)), std(out_m.dff_rew_willSwitch{ne}(logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_rew_willSwitch{ne}(logical(out_m.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.maleC})
    title('CS plus preferring')
    xlabel('Time from reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out_f.dff_rew_switch{ne},2));
    subplot(2,2,2); hold on, box off
    lh{1}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willStay{ne}(logical(out_f.negKernel{ne}),:)), std(out_f.dff_nrew_willStay{ne}(logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_nrew_willStay{ne}(logical(out_f.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.femaleC});
    lh{2}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willSwitch{ne}(logical(out_f.negKernel{ne}),:)), std(out_f.dff_nrew_willSwitch{ne}(logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_nrew_willSwitch{ne}(logical(out_f.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.femaleC});
     lh{3}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willStay{ne}(logical(out_m.negKernel{ne}),:)), std(out_m.dff_nrew_willStay{ne}(logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_nrew_willStay{ne}(logical(out_m.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.maleC});
    lh{4}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willSwitch{ne}(logical(out_m.negKernel{ne}),:)), std(out_m.dff_nrew_willSwitch{ne}(logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_nrew_willSwitch{ne}(logical(out_m.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.maleC});
    title('CS plus preferring')
    xlabel('Time from no reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    subplot(2,2,3); hold on, box off
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willStay{ne}(~logical(out_f.negKernel{ne}),:)), std(out_f.dff_rew_willStay{ne}(~logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_rew_willStay{ne}(~logical(out_f.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_f.dff_rew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:)), std(out_f.dff_rew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_rew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.femaleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willStay{ne}(~logical(out_m.negKernel{ne}),:)), std(out_m.dff_rew_willStay{ne}(~logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_rew_willStay{ne}(~logical(out_m.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.maleC})
    shadedErrorBar(xaxis,mean(out_m.dff_rew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:)), std(out_m.dff_rew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_rew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.maleC})
  title('CS minus preferring')
    xlabel('Time from reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out_f.dff_rew_switch{ne},2));
    subplot(2,2,4); hold on, box off
    lh{1}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willStay{ne}(~logical(out_f.negKernel{ne}),:)), std(out_f.dff_nrew_willStay{ne}(~logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_nrew_willStay{ne}(~logical(out_f.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.femaleC});
    lh{2}=shadedErrorBar(xaxis,mean(out_f.dff_nrew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:)), std(out_f.dff_nrew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:))./sqrt(size(out_f.dff_nrew_willSwitch{ne}(~logical(out_f.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.femaleC});
    lh{3}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willStay{ne}(~logical(out_m.negKernel{ne}),:)), std(out_m.dff_nrew_willStay{ne}(~logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_nrew_willStay{ne}(~logical(out_m.negKernel{ne}),:),1)),'LineProps',{'Color',plotParams.maleC});
    lh{4}=shadedErrorBar(xaxis,mean(out_m.dff_nrew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:)), std(out_m.dff_nrew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:))./sqrt(size(out_m.dff_nrew_willSwitch{ne}(~logical(out_m.negKernel{ne}),:),1)),'LineProps',{'--','Color',plotParams.maleC});
    title('CS minus preferring')
    xlabel('Time from no reward cue (s)')
    ylabel('Mean fluorescence (Z-score)')
    legend([lh{1}.mainLine lh{2}.mainLine lh{3}.mainLine lh{4}.mainLine],{'female: Stay', 'female:Switch','male: Stay', 'male:Switch'})
    
    % Plot mean activity across time for each cell 
    figure()
    subplot(1,2,1);hold on, box off
    scatter(mean(out_f.dff_rew_willStay{ne},2), mean(out_f.dff_rew_willSwitch{ne},2),12,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5) 
    scatter(mean(out_m.dff_rew_willStay{ne},2), mean(out_m.dff_rew_willSwitch{ne},2),12,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5) 

    plot([-.5 1],[-.5 1],':','Color',[.7 .7 .7])
    xlabel('Average rew -> stay')
    ylabel('Average rew -> switch')
    axis square
    subplot(1,2,2);hold on, box off
    plot([-.5 1],[-.5 1],':','Color',[.7 .7 .7])
    scatter(mean(out_f.dff_nrew_willStay{ne},2), mean(out_f.dff_nrew_willSwitch{ne},2),12,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5) 
    scatter(mean(out_m.dff_nrew_willStay{ne},2), mean(out_m.dff_nrew_willSwitch{ne},2),12,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5) 

    xlabel('Average no rew -> stay')
    ylabel('Average no rew -> switch')
    axis square
    
    gcamp = cat(1,mean(out_f.dff_rew_willStay{ne},2),mean(out_f.dff_rew_willSwitch{ne},2),mean(out_f.dff_nrew_willStay{ne},2),mean(out_f.dff_nrew_willSwitch{ne},2),...
        mean(out_m.dff_rew_willStay{ne},2),mean(out_m.dff_rew_willSwitch{ne},2),mean(out_m.dff_nrew_willStay{ne},2),mean(out_m.dff_nrew_willSwitch{ne},2));
    X = table(gcamp);
    X.Sex = cat(1,ones(size(out_f.dff_rew_willStay{ne},1),1),ones(size(out_f.dff_rew_willSwitch{ne},1),1),ones(size(out_f.dff_nrew_willStay{ne},1),1),ones(size(out_f.dff_nrew_willSwitch{ne},1),1),...
        zeros(size(out_m.dff_rew_willStay{ne},1),1),zeros(size(out_m.dff_rew_willSwitch{ne},1),1),zeros(size(out_m.dff_nrew_willStay{ne},1),1),zeros(size(out_m.dff_nrew_willSwitch{ne},1),1));
    X.Stay = cat(1,ones(size(out_f.dff_rew_willStay{ne},1),1),zeros(size(out_f.dff_rew_willSwitch{ne},1),1),ones(size(out_f.dff_nrew_willStay{ne},1),1),zeros(size(out_f.dff_nrew_willSwitch{ne},1),1),...
        ones(size(out_m.dff_rew_willStay{ne},1),1),zeros(size(out_m.dff_rew_willSwitch{ne},1),1),ones(size(out_m.dff_nrew_willStay{ne},1),1),zeros(size(out_m.dff_nrew_willSwitch{ne},1),1));
    X.Reward = cat(1,ones(size(out_f.dff_rew_willStay{ne},1),1),ones(size(out_f.dff_rew_willSwitch{ne},1),1),zeros(size(out_f.dff_nrew_willStay{ne},1),1),zeros(size(out_f.dff_nrew_willSwitch{ne},1),1),...
        ones(size(out_m.dff_rew_willStay{ne},1),1),ones(size(out_m.dff_rew_willSwitch{ne},1),1),zeros(size(out_m.dff_nrew_willStay{ne},1),1),zeros(size(out_m.dff_nrew_willSwitch{ne},1),1));

    X.Sex = nominal(X.Sex);
    X.Stay = nominal(X.Stay);
    X.Reward = nominal(X.Reward);
    glme = fitglme(X,'gcamp~Sex*Stay*Reward');
    
    
    
    % Plot max activity across time for each cell 
    figure()
    subplot(1,2,1);hold on, box off
    [~,peakidx] = max(abs(out_f.dff_rew_willStay{ne}),[],2);
    [~,peakidx_2] = max(abs(out_f.dff_rew_willSwitch{ne}),[],2);
    max_rew_stay = arrayfun(@(x,y) out_f.dff_rew_willStay{ne}(x,y),[1:numel(peakidx)]',peakidx); 
    max_rew_switch = arrayfun(@(x,y) out_f.dff_rew_willSwitch{ne}(x,y),[1:numel(peakidx_2)]',peakidx_2);
    
    scatter(max_rew_stay,max_rew_switch,12,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.7) 
    
    [~,peakidx] = max(abs(out_m.dff_rew_willStay{ne}),[],2);
    [~,peakidx_2] = max(abs(out_m.dff_rew_willSwitch{ne}),[],2);
    max_rew_stay = arrayfun(@(x,y) out_m.dff_rew_willStay{ne}(x,y),[1:numel(peakidx)]',peakidx); 
    max_rew_switch = arrayfun(@(x,y) out_m.dff_rew_willSwitch{ne}(x,y),[1:numel(peakidx_2)]',peakidx_2);
    
    scatter(max_rew_stay,max_rew_switch,12,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.7) 
    
    plot([-.5 1],[-.5 1],':','Color',[.7 .7 .7],'LineWidth',1)

    xlabel('Peak rew -> stay')
    ylabel('Peak rew -> switch')
    axis square
    
     subplot(1,2,2);hold on, box off
    [~,peakidx] = max(abs(out_f.dff_nrew_willStay{ne}),[],2);
    [~,peakidx_2] = max(abs(out_f.dff_nrew_willSwitch{ne}),[],2);
    max_nrew_stay = arrayfun(@(x,y) out_f.dff_nrew_willStay{ne}(x,y),[1:numel(peakidx)]',peakidx); 
    max_nrew_switch = arrayfun(@(x,y) out_f.dff_nrew_willSwitch{ne}(x,y),[1:numel(peakidx_2)]',peakidx_2);
    
    scatter(max_nrew_stay,max_nrew_switch,12,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.7) 
    
    [~,peakidx] = max(abs(out_m.dff_nrew_willStay{ne}),[],2);
    [~,peakidx_2] = max(abs(out_m.dff_nrew_willSwitch{ne}),[],2);
    max_nrew_stay = arrayfun(@(x,y) out_m.dff_nrew_willStay{ne}(x,y),[1:numel(peakidx)]',peakidx); 
    max_nrew_switch = arrayfun(@(x,y) out_m.dff_nrew_willSwitch{ne}(x,y),[1:numel(peakidx_2)]',peakidx_2);
    
    scatter(max_nrew_stay,max_nrew_switch,12,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.7) 
    
    plot([-.5 1],[-.5 1],':','Color',[.7 .7 .7],'LineWidth',1)

    xlabel('Peak nrew -> stay')
    ylabel('Peak nrew -> switch')
    axis square
    
    %% Plot histogram of correlation coefficients 
    figure()
    subplot(2,2,1); hold on
    [a,b] = histcounts(out_f.corr_rew_willStay{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.femaleC)
    [a,b] = histcounts(out_m.corr_rew_willStay{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.maleC)
    ylabel('Proportion of neurons')
    xlabel('Correlation coefficient')
    title('Rew -> Stay X latency')
     subplot(2,2,2); hold on
    [a,b] = histcounts(out_f.corr_rew_willSwitch{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.femaleC)
    [a,b] = histcounts(out_m.corr_rew_willSwitch{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.maleC)
    ylabel('Proportion of neurons')
    xlabel('Correlation coefficient')
    title('Rew -> Switch X latency')
    
    subplot(2,2,3); hold on
    [a,b] = histcounts(out_f.corr_nrew_willStay{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.femaleC)
    [a,b] = histcounts(out_m.corr_nrew_willStay{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.maleC)
    ylabel('Proportion of neurons')
    xlabel('Correlation coefficient')
    title('No Rew -> Stay X latency')
     subplot(2,2,4); hold on
    [a,b] = histcounts(out_f.corr_nrew_willSwitch{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.femaleC)
    [a,b] = histcounts(out_m.corr_nrew_willSwitch{ne},[-1:.1:1],'normalization','probability');
    stairs(b(1:end-1),a,'Color',plotParams.maleC)
    ylabel('Proportion of neurons')
    xlabel('Correlation coefficient')
    title('No Rew -> Switch X latency')
end
end

function out=extractStaySwitch(recs,rawFlag,frameRate,zscoreFlag,sigVer,sigEvent,qFile,events,eventDur,kernelVer)

findSig = 0;
fbasename = fullfile(whereAreWe('imaging'));
nCVFolds = 5;

% Load basis set
fbasename_bs = fullfile(whereAreWe('figurecode'),'general_code', 'basis_sets');
[~, ~, ~, ~, ~, ~, bsIDs,~,~] = getEvents(sigVer,frameRate);

load(fullfile(fbasename_bs, ['bs_' bsIDs '.mat']))
bs = (full(eval(['bs_' bsIDs])));



a = 0.01; 


% Initialize summary variables
out.dff_rew_stay = cell(numel(events));
out.dff_nrew_stay = cell(numel(events));
out.dff_rew_willStay = cell(numel(events));
out.dff_nrew_willStay = cell(numel(events));


out.dff_rew_switch = cell(numel(events));
out.dff_nrew_switch = cell(numel(events));
out.dff_rew_willSwitch = cell(numel(events));
out.dff_nrew_willSwitch = cell(numel(events));


out.corr_rew_stay = cell(numel(events));
out.corr_nrew_stay = cell(numel(events));
out.corr_rew_willStay = cell(numel(events));
out.corr_nrew_willStay = cell(numel(events));


out.corr_rew_switch = cell(numel(events));
out.corr_nrew_switch = cell(numel(events));
out.corr_rew_willSwitch = cell(numel(events));
out.corr_nrew_willSwitch = cell(numel(events));

out.negKernel = cell(numel(events));

%% Extract event-triggered activity for neurons significant in sigVer regression for sigEvent
roiIdx = 1;
for na = 1:numel(recs)
    clear pmat
    
    %% Load results for regression to select significantly encoding neurons for all events
    load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,sigVer)),'pvals','b'); %load regression results
    load(fullfile(fbasename,recs{na},['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat'])); % load events
    %% Load value and dff
%     try load(fullfile(fbasename,recs{na},sprintf('valueData_fs%d_raw%d_nfolds%d_%s.mat',frameRate,rawFlag,nCVFolds,qFile)));
%     catch
%         fprintf('extracting value for %s \n',recs{na});
%         pause(.01)
%         setupValueModel(recs{na}, rawFlag, qFile, frameRate, tdtEvents, nCVFolds);
%         load(fullfile(fbasename, recs{na}, sprintf('valueData_fs%s_raw%s_nfolds%s_%s', num2str(frameRate), num2str(rawFlag), num2str(nCVFolds), qFile)));
%     end
    
    % Find significant neurons
    try load(fullfile(fbasename,recs{na},'activeNeurons.mat'))
    catch
        activeIdx = findActiveNeurons(recs{na});
    end
    
    pvals = cell2mat(pvals');
    if ~findSig
        % all neurons
        sigIdx = activeIdx;
    else
        % significant neurons from regression
        try
            sigIdx = find(pvals(activeIdx,sigEvent)<a);
        catch
            keyboard
        end
    end
    
    % load dff
    load(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'),sprintf('dff_%d_raw%d',frameRate,rawFlag));
    dff = eval(sprintf('dff_%d_raw%d',frameRate,rawFlag));
    eval(sprintf('clear dff_%d_raw%d',frameRate,rawFlag));
    
    if zscoreFlag
        dff = nanzscore(dff,1);
    end
    
    reward = tdtEvents.reward;
    stay = [NaN tdtEvents.choice(1:end-1)==tdtEvents.choice(2:end)];
    willStay = [stay(2:end) NaN];
    latency = tdtEvents.nosePokeEntryTimestamps-tdtEvents.trialStartTimestamps;
    latency = [latency(2:end),NaN]';
    for ne = 1:numel(events)
       if contains(events{ne},'outcome')
           reward = tdtEvents.reward;
       else
           reward = [NaN tdtEvents.reward(1:end-1)]; %switch to previous reward if pre-outcome event
       end
       
        
       thisEvent = eval(sprintf('tdtEvents.%s',events{ne}));
       for nn = sigIdx'
           load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,kernelVer)),'b','con_iden'); %load regression results for kernel
           load(fullfile(fbasename, recs{na}, sprintf('M_linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,kernelVer)),'x_basic'); %load regression results
           intercept = b{nn}(1);
           bb = b{nn}(2:end);
           x_basic = nanzscore(x_basic,1);
           xidx = (ismember(con_iden,1:numel(unique(con_iden))-2));
           y_pred = intercept+x_basic(:,xidx)*bb(xidx); % predicted values
           try
           y_resid = dff(:,nn)-y_pred;
           catch
               keyboard
           end
           
           load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,sigVer)),'b'); %load regression results for significance 
           
           out.latency{roiIdx} = latency;
           out.willStay{roiIdx} = willStay;
           out.reward{roiIdx} = tdtEvents.reward;
           out.choice{roiIdx} = [tdtEvents.ipsiChoice(2:end),NaN];
           out.trialStart{roiIdx} = tdtEvents.trialStartTimestamps; 
           out.nosePokeEntry{roiIdx} = tdtEvents.nosePokeEntryTimestamps;
           out.outcomeTime{roiIdx} = zeros(size(tdtEvents.outcome));
           out.outcomeTime{roiIdx}(tdtEvents.reward==1) = tdtEvents.CSPlusTimestamps;
           out.outcomeTime{roiIdx}(tdtEvents.reward==0&tdtEvents.choice~=-1) = tdtEvents.CSMinusTimestamps;

           clear thistrials thistrials_resid
           % Extract activity around event
           for nt = 1:numel(thisEvent)
               if isnan(thisEvent(nt))
                   thistrials(nt,:) = nan(size(-eventDur(ne,1)*frameRate:eventDur(ne,2)*frameRate));
                   thistrials_resid(nt,:) = nan(size(-eventDur(ne,1)*frameRate:eventDur(ne,2)*frameRate));
               else
                   try
                       thistrials(nt,:) = dff(thisEvent(nt)-eventDur(ne,1)*frameRate:thisEvent(nt)+eventDur(ne,2)*frameRate,nn);
                       thistrials_resid(nt,:) = y_resid(thisEvent(nt)-eventDur(ne,:)*frameRate:thisEvent(nt)+eventDur(ne,2)*frameRate);
                   catch
                       thistrials(nt,:) = cat(1,dff(thisEvent(nt)-eventDur(ne,1)*frameRate:end,nn),zeros(thisEvent(nt)+eventDur(ne,2)*frameRate-size(dff,1),1));
                       thistrials_resid(nt,:) = cat(1,y_resid(thisEvent(nt)-eventDur(ne,1)*frameRate:end),zeros(thisEvent(nt)+eventDur(ne,2)*frameRate-size(dffi,1),1)); 
                   end
               end
           end
           out.thistrials_resid{roiIdx} = thistrials_resid;
           out.thistrials_orig{roiIdx} = thistrials;
           if findSig
               if contains(events{ne},'outcome')
                   p = ranksum(nanmean(thistrials(willStay==1,:),2),nanmean(thistrials(willStay==0,:),2));
               else
                   p = ranksum(nanmean(thistrials(stay==1,:),2),nanmean(thistrials(stay==0,:),2));
               end
               if p<a
               else
                   continue
               end
           end
           
           
           if size(thistrials,1)==numel(tdtEvents.choice)
               % Extract event activity for each reward sequence
               out.dff_rew_stay{ne} = cat(1,out.dff_rew_stay{ne},nanmean(thistrials(reward==1&stay==1,:)));
               out.dff_nrew_stay{ne} = cat(1,out.dff_nrew_stay{ne},nanmean(thistrials(reward==0&stay==1,:)));
               out.dff_rew_switch{ne} = cat(1,out.dff_rew_switch{ne},nanmean(thistrials(reward==1&stay==0,:)));
               out.dff_nrew_switch{ne} = cat(1,out.dff_nrew_switch{ne},nanmean(thistrials(reward==0&stay==0,:)));
               
               out.dff_rew_willStay{ne} = cat(1,out.dff_rew_willStay{ne},nanmean(thistrials(reward==1&willStay==1,:),1));
               out.dff_nrew_willStay{ne} = cat(1,out.dff_nrew_willStay{ne},nanmean(thistrials(reward==0&willStay==1,:),1));
               out.dff_rew_willSwitch{ne} = cat(1,out.dff_rew_willSwitch{ne},nanmean(thistrials(reward==1&willStay==0,:),1));
               out.dff_nrew_willSwitch{ne} = cat(1,out.dff_nrew_willSwitch{ne},nanmean(thistrials(reward==0&willStay==0,:),1));
            
               out.dff_rew_stay_trials{ne}{roiIdx} = thistrials(reward==1&stay==1,:);
               out.dff_nrew_stay_trials{ne}{roiIdx} = thistrials(reward==0&stay==1,:);
               out.dff_rew_switch_trials{ne}{roiIdx} = thistrials(reward==1&stay==0,:);
               out.dff_nrew_switch_trials{ne}{roiIdx} = thistrials(reward==0&stay==0,:);
               
               out.dff_rew_willStay_trials{ne}{roiIdx} = thistrials(reward==1&willStay==1,:);
               out.dff_nrew_willStay_trials{ne}{roiIdx} = thistrials(reward==0&willStay==1,:);
               out.dff_rew_willSwitch_trials{ne}{roiIdx} = thistrials(reward==1&willStay==0,:);
               out.dff_nrew_willSwitch_trials{ne}{roiIdx} = thistrials(reward==0&willStay==0,:);
               
               % Correlation between event activity and latency
               out.corr_rew_stay{ne} = cat(1,out.corr_rew_stay{ne},corr(latency(reward==1&stay==1),nanmean(thistrials(reward==1&stay==1,:),2),'rows','complete'));
               out.corr_nrew_stay{ne} = cat(1,out.corr_nrew_stay{ne},corr(latency(reward==0&stay==1),nanmean(thistrials(reward==0&stay==1,:),2),'rows','complete'));
               
               out.corr_rew_switch{ne} = cat(1,out.corr_rew_switch{ne},corr(latency(reward==1&stay==0),nanmean(thistrials(reward==1&stay==0,:),2),'rows','complete'));
               out.corr_nrew_switch{ne} = cat(1,out.corr_nrew_switch{ne},corr(latency(reward==0&stay==0),nanmean(thistrials(reward==0&stay==0,:),2),'rows','complete'));
               
               out.corr_rew_willStay{ne} = cat(1,out.corr_rew_willStay{ne},corr(latency(reward==1&willStay==1),nanmean(thistrials(reward==1&willStay==1,:),2),'rows','complete'));
               out.corr_nrew_willStay{ne} = cat(1,out.corr_nrew_willStay{ne},corr(latency(reward==0&willStay==1),nanmean(thistrials(reward==0&willStay==1,:),2),'rows','complete'));
               
               out.corr_rew_willSwitch{ne} = cat(1,out.corr_rew_willSwitch{ne},corr(latency(reward==1&willStay==0),nanmean(thistrials(reward==1&willStay==0,:),2),'rows','complete'));
               out.corr_nrew_willSwitch{ne} = cat(1,out.corr_nrew_willSwitch{ne},corr(latency(reward==0&willStay==0),nanmean(thistrials(reward==0&willStay==0,:),2),'rows','complete'));

                    
           else
               keyboard
           end
           
           % Positive or negatively encoding sigEvent?
           if ~iscell(bsIDs)
               thisBeta = b{nn}(2+(sigEvent-1)*size(bs,2):sigEvent*size(bs,2)+1);
           else
               keyboard
           end
           thisKernel = bs*thisBeta;
           if trapz(thisKernel)<0
               out.negKernel{ne} = cat(1,out.negKernel{ne},1);
           elseif trapz(thisKernel)>0
               out.negKernel{ne} = cat(1,out.negKernel{ne},0);
           else
               keyboard
           end
           
           
           
           

           roiIdx=roiIdx+1;

       end
    end  
end  
    
end



