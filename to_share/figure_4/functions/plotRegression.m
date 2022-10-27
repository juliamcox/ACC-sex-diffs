function plotRegression(recs, sigVer, plotVer, frameRate, rawFlag, zscoreFlag, savename, sigEvents,plotEvents,norm,histLen,sigLevel,aoc_event,cohort)

%aids: cell array of animal IDs
%recs: cell array of recording names
%sigVer: version of regression to use for significance (Basic3: all events;fig3_outcome;fig3_choice)
%plotVer: version of regression to plot (fHist to plot gcamp)
%frameRate
%rawFlag: were the fluorescence traces raw or denoised?
%zscoreFlag: were the fluoresence traces zscored
%savename: file name for saving 
%norm: none, zscore, peak
%histLen: needs to be entered only if plotting fHist 

% sigEvents: which events to use to determine significance 
% plotEvents: which epoch to plot; if whicEpoch = 'fHist', plotEvents should be cell array of event names 
% when plotEvents is a cell array of event names, histLength is a number of epochs X 2 matrix with plot start and end times for each event (in seconds)

if iscell(plotEvents) && nargin < 11
    error('Need to input plot lengths for each epoch in plotEvents')
end
if nargin < 14
    cohort = 'ACC_DMS_imaging';
end

numPlots = 8;
plotRaw = 1;

counter = 0;
counter2=0;
fbasename = fullfile(whereAreWe('imaging'),cohort);
fbasename_bs = fullfile(whereAreWe('bucket'),'DMS_Bandit', 'basis sets');
a = sigLevel;
if contains(plotVer,'fHist')
    cons = 1:numel(plotEvents);
else
    cons = getEvents(plotVer,frameRate);

end

% initialize variables
pmat_all = cell(size(recs));
neuronIdx = cell(size(recs));
b_all = cell(size(recs));
b_bs_all = cell(size(cons));

if exist('histLen')
    histLenF = histLen.*frameRate;
end

for na = 1:numel(recs)
    clear pmat
    
    % load regression results
    %cd(fullfile(fbasename, recs{na}));
    info = matfile(fullfile(fbasename,recs{na},'info.mat'));
    frameRateO = info.frameRate;
    droppedIdx = info.droppedIdx;
    nframesLog = info.nframesLog;
    
    %% Load results for regression to select significantly encoding neurons for all events
    try
        load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,sigVer)),'pvals','con_iden'); %load regression results
    catch
        load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig4_%s.mat', frameRate,rawFlag,zscoreFlag,sigVer)),'pvals','con_iden'); %load regression results
    end
    % number of events
    numEvents = numel(sigEvents);
    clear pmat
    for nn = 1:numel(pvals)
        pmat(nn,1:size(pvals{nn},2)) = pvals{nn} < a;
    end
    [r,c] =(find(pmat(:,sigEvents)));
    idx = unique(r);
    
    try load(fullfile(fbasename,recs{na},'activeNeurons.mat'))
    catch
        activeIdx = findActiveNeurons(recs{na},fbasename);
    end
    
    pmat_all{na} = pmat;
    neuronIdx{na} = idx(ismember(idx,activeIdx));
    %neuronIdx{na} = idx;
    %save(fullfile(fbasename, recs{na}, sprintf('pmat_linReg_fig2_%s_fs%d_raw%d_zscore%d_%s', sigVer, frameRate, rawFlag, zscoreFlag)), 'pmat')
    counter = counter+length(activeIdx);
    counter2 = counter2+length(pvals); 
    
    
    %% Load results for plotVer and extract predictors
    
    if contains(plotVer, 'fHist')
        % Load event times
        load(fullfile(fbasename,recs{na},['behavEventsRegression_bs', num2str(1000/round(frameRate)) '.mat']));
        % Load gcamp
        
        try
            load(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'),sprintf('dff_%d_raw%d',frameRate,plotRaw));
            dff = eval(sprintf('dff_%d_raw%d',frameRate,plotRaw));
            eval(sprintf('clear dff_%d_raw%d',frameRate,plotRaw));
         catch
            load(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'),'neuron');
            if plotRaw %if rawFlag use the raw fluorescence
                dff = neuron.C_raw;
            else %otherwise use denoised fluorescence
                dff = neuron.C;
            end
            
            if size(dff,2) > size(dff,1)
                dff = dff';
            end
            %if droppped frames had been removed for cnmfe, add back in NaN
            if droppedIdx == 0
                droppedIdx = [];
            end
            if numel(droppedIdx)>0 numel(droppedIdx)+nframesLog > length(dff)
                temp = zeros(numel(droppedIdx)+nframesLog,size(dff,2));
                temp(droppedIdx,:) = NaN;
                temp(setdiff(1:numel(droppedIdx)+nframesLog,droppedIdx),:) = dff;
                dff = temp;
            elseif numel(droppedIdx) == 0 && nframesLog ~= length(dff)
                keyboard %this is probably T55 7/24 where I had to cut two frames to get inscopix to export tiff
            end
            %downsample if necessary
            if frameRate~=frameRateO
                fprintf('Downsampling....\n')
                %dff = downsampleImaging(frameRateO, frameRate,dff);
                dff = resample(dff,round(frameRate),round(frameRateO));
            end
            eval(sprintf('dff_%d_raw%d = dff;', frameRate,plotRaw));
            save(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'), sprintf('dff_%d_raw%d', frameRate,plotRaw), '-append');
        end
        
        if size(dff,1)<size(dff,2)
            dff = dff';
        end
        dff = nanzscore(dff);
        
        if ~iscell(plotEvents)
            error('plotEvents must contain event names to plot gcamp histograms')
        end
        for ne = 1:numel(plotEvents)
            if contains(plotEvents{ne},'ipsiLeverPress')
                if info.ipsiSide == 1
                    thisEvent = tdtEvents.LLeverPress;
                else
                    thisEvent = tdtEvents.RLeverPress;
                end
            elseif contains(plotEvents{ne}, 'contraLeverPress')
                if info.ipsiSide == 1
                    thisEvent = tdtEvents.RLeverPress;
                else
                    thisEvent = tdtEvents.LLeverPress;
                end
            else
                thisEvent = eval(sprintf('tdtEvents.%s',plotEvents{ne}));
            end
            
            for nr = neuronIdx{na}'
                for nt = 1:numel(thisEvent)
                    thisidx = thisEvent(nt)+histLenF(ne,1):thisEvent(nt)+histLenF(ne,2);
                    try
                        if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                            thistrial = [dff(thisEvent(nt)+histLenF(1):end,nr) zeros(thisidx(end)-length(dff),1)];
                        else thistrial = dff(thisidx,nr);
                        end
                        try
                            thishist(nt,:) = thistrial;
                        catch
                            keyboard
                        end
                    catch
                        keyboard
                    end
                end
                b_bs_all{ne} = cat(1,b_bs_all{ne},nanmean(thishist));
            end
            
        end
        
        
    else
        
%        load(fullfile(fbasename, 'predictors', [plotVer '.mat']), 'bsIDs','cons','time_back_orig','time_forward_orig','con_shift') % Load regression parameters for plotVer
        [cons, con_shift, time_back_orig, time_forward_orig, ~, ~, bsIDs,~,~] = getEvents(plotVer,frameRate);

        try
            load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s.mat', frameRate,rawFlag,zscoreFlag,plotVer)),'con_iden','b'); % Load regression coefficients for plotVer
        catch
            keyboard
            load(fullfile(fbasename, recs{na}, sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig4_%s.mat', frameRate,rawFlag,zscoreFlag,plotVer)),'con_iden','b'); % Load regression coefficients for plotVer
        end
        
        
        % load basis set
        if iscell(bsIDs)
            for nb = 1:numel(bsIDs)
                load(fullfile(fbasename_bs, ['bs_' bsIDs{nb} '.mat']))
                bs{nb} = (full(eval(['bs_' bsIDs{nb}])));
            end
        else
            load(fullfile(fbasename_bs, ['bs_' bsIDs '.mat']))
            bs = (full(eval(['bs_' bsIDs])));
        end
        b = cellfun(@(x) x(2:end),b,'UniformOutput',false);
        b_all{na} = cell2mat(b(neuronIdx{na}));
        
        % calculate kernel: extract coefficients for each event and multiply with basis set and sum
        for nn = 1:numel(neuronIdx{na})
            for ne = plotEvents
                thisWeights = b{neuronIdx{na}(nn)}(con_iden==ne);
                if iscell(bs)
                    tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*bs{ne},2)';
                else
                    tempWeights = sum(repmat(thisWeights',size(bs,1),1).*bs,2)';
                end
                b_bs_all{ne} = cat(1,b_bs_all{ne}, tempWeights);
                clear tempWeights
            end
        end
        
    end
    
end

%% Plot


fprintf('%d \n',size(cell2mat(b_all'),2)/counter)
fprintf('%d \n', size(cell2mat(b_all'),2))
fprintf('%d \n', counter)
fprintf('%d \n', counter2)

%% Plot based on area under the curve
if nargin>12
    
    fprintf('load sort \n')
    keyboard
    negIdx = find(aoc_event<0);
    posIdx = find(aoc_event>0);
    
    if contains(sigVer,'outcome')
        plot1 = cell2mat(b_bs_all(2));
    else
        plot1 = cell2mat(b_bs_all(1));
    end
    plot1 = plot1(posIdx,:);
    
    [~,maxI] = max((plot1),[],2);
    [temp, maxI] = sort(maxI);
    plot1 = cell2mat(b_bs_all);
    plot1 = plot1(posIdx,:);
    try
        plot1 = plot1(maxI_pos,:);
    catch
        plot1 = plot1(maxI,:);
    end
    
    if contains(sigVer,'outcome')
        plot2 = cell2mat(b_bs_all(1));
    else
        plot2 = cell2mat(b_bs_all(2));
    end
    plot2 = plot2(negIdx,:);
     
    [~,maxI] = max((plot2),[],2);
   [temp, maxI] = sort(maxI);
    plot2 = cell2mat(b_bs_all);
    plot2 = plot2(negIdx,:);
    try
        plot2 = plot2(maxI_neg,:);
        fprintf('No sort loaded \n');
    catch
        plot2 = plot2(maxI,:);
    end
    if contains(sigVer,'outcome')
    plotAll = cat(1,plot2,plot1);
    else
        plotAll = cat(1,plot1,plot2);
    end
        
    
    if ~contains(plotVer,'fHist')
        switch norm
            case 'none'
                maxC = max(max(plotAll));
                minC = min(min(plotAll));
                 minC = prctile(max(cell2mat(b_bs_all)),10);
            maxC = prctile(max(cell2mat(b_bs_all)),90); 
                if abs(minC)>abs(maxC)
                    maxC = -minC;
                else
                    minC = -maxC;
                end
                temp2 = plotAll;
            case 'zscore'
                temp2 = plotAll;
                mu = nanmean(temp2,2);
                sigma = nanstd(temp2,[],2);
                temp2 = (temp2-repmat(mu,1,size(temp2,2)))./repmat(sigma,1,size(temp2,2));
                maxC = 6;
                minC = -6;
            case 'peak'
                temp2 = plotAll;
                temp2=temp2./repmat(max((temp2),[],2),1,size(temp,2));
                maxC = 1;
                minC = -1;
        end
    else
        maxC = max(max(max(plotAll)));
        minC = min(min(min(plotAll)));
          minC = prctile(max(cell2mat(b_bs_all)),5);
            maxC = prctile(max(cell2mat(b_bs_all)),95); 
        if abs(minC)>abs(maxC)
            maxC = -minC;
        else
            minC = -maxC;
        end
        temp2 = (plotAll);
    end
    
    
    minC = -3;
    maxC = 3;
    
    screensize = get( groot, 'Screensize' );
    screensize(4) = screensize(4)/1.5;
    screensize(3) = screensize(3)/1.5;
    %screensize(2) = screensize(2)+screensize(4)-100;
    figure('Position', screensize)
    
    if ~contains(plotVer,'fHist')
        x=1;
        for ne = plotEvents
            idx = x:size(b_bs_all{ne},2)+x-1;
            a=subplot(1,numPlots,ne);
            if numel(time_forward_orig) > 1
                xaxis = linspace(-time_back_orig(ne),time_forward_orig(ne), size(b_bs_all{ne},2));
            else
                if con_shift(ne) == 0
                    xaxis = linspace(-2,6,size(b_bs_all{ne},2));
                else
                    xaxis = linspace(0,8,size(b_bs_all{ne},2));
                end
            end
            imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(:,idx), [minC maxC]);
            colormap(flipud(red2blue));
            %      colormap(bluewhitered);
            
            box off
            title(cons{ne});
            if ne == plotEvents(end)
                b=smallcolorbar(gca);
                b.Ticks = [minC 0 maxC];
                b.YLabel.String = 'Kernel value';
            end
            if ne ~= plotEvents(1)
                a.YTick = [];
            else
                ylabel('Neuron number (sorted by time of peak)')
            end
            x = x+size(b_bs_all{ne},2);
        end
        
    else
        x=1;
        for ne = 1:numel(plotEvents)
            idx = x:size(b_bs_all{ne},2)+x-1;
            a=subplot(1,numPlots,ne);
            xaxis = linspace(histLen(ne,1),histLen(ne,2),size(b_bs_all{ne},2));
            
            imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(:,idx),[minC maxC]);
            colormap(flipud(red2blue));
            %colormap(bluewhitered);
            box off
            title(plotEvents{ne});
            if ne == numel(plotEvents)
                b=smallcolorbar(gca);
                b.Ticks = [minC 0 maxC];
                b.YLabel.String = 'Kernel value';
            end
            if ne ~= 1
                a.YTick = [];
            else
                ylabel('Neuron number (sorted by time of peak)')
            end
            x = x+size(b_bs_all{ne},2);
        end
    end
    if contains(sigVer,'outcome')
    hold on, plot(xaxis,ones(size(xaxis)).*numel(negIdx),'k')
    else
            hold on, plot(xaxis,ones(size(xaxis)).*numel(posIdx),'k')
    end
end

%% Plot heatmap of kernels for all neurons sorted by time of peak
[~,maxI_orig] = max(cell2mat(b_bs_all),[],2);
[~, maxI] = sort(maxI_orig);

% if contains(plotVer,'fHist')
%     load(fullfile(whereAreWe('bucket'),'DMS_Bandit','Basic3_maxI'));
%     keyboard
% end

%normalize regression coefficients 
if ~contains(plotVer,'fHist')
    switch norm
        case 'none'
            temp2 = cell2mat(b_bs_all);
            maxC = max(max(cell2mat(b_bs_all)));
            minC = min(min(cell2mat(b_bs_all)));
            minC = prctile(max(cell2mat(b_bs_all)),5);
            maxC = prctile(max(cell2mat(b_bs_all)),95); 
            if abs(minC)>abs(maxC)
                maxC = -minC;
            else
                minC = -maxC;
            end
        case 'zscore'
            temp2 = cell2mat(b_bs_all);
            mu = nanmean(temp2,2);
            sigma = nanstd(temp2,[],2);
            temp2 = (temp2-repmat(mu,1,size(temp2,2)))./repmat(sigma,1,size(temp2,2));
            maxC = 6;
            minC = -6;
        case 'peak'
            temp2 = cell2mat(b_bs_all);
            temp2=temp2./repmat(max((temp2),[],2),1,size(temp,2));
            maxC = 1;
            minC = -1;
    end
else
    maxC = max(max(max(cell2mat(b_bs_all))));
    minC = min(min(min(cell2mat(b_bs_all))));
    if abs(minC)>abs(maxC)
        maxC = -minC;
    else
        minC = -maxC;
    end
    temp2 = cell2mat(b_bs_all); 
end




screensize = get( groot, 'Screensize' );
screensize(4) = screensize(4)/1.5;
screensize(3) = screensize(3)/1.5;
%screensize(2) = screensize(2)+screensize(4)-100;
figure('Position', screensize)

if ~contains(plotVer,'fHist')
    x=1;
    for ne = plotEvents
        idx = x:size(b_bs_all{ne},2)+x-1;
        a=subplot(1,numPlots,ne);
        if numel(time_forward_orig) > 1
            xaxis = linspace(-time_back_orig(ne),time_forward_orig(ne), size(b_bs_all{ne},2));
        else
            if con_shift(ne) == 0
                xaxis = linspace(-2,6,size(b_bs_all{ne},2));
            else
                xaxis = linspace(0,8,size(b_bs_all{ne},2));
            end
        end
        imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(maxI,idx), [minC maxC]);
        colormap(flipud(red2blue));
          %      colormap(bluewhitered);

        box off
        title(cons{ne});
        if ne == plotEvents(end)
            b=smallcolorbar(gca);
            b.Ticks = [minC 0 maxC];
            b.YLabel.String = 'Kernel value';
        end
        if ne ~= plotEvents(1)
            a.YTick = [];
        else
            ylabel('Neuron number (sorted by time of peak)')
        end
        x = x+size(b_bs_all{ne},2);
    end
    
else
    x=1;
    for ne = 1:numel(plotEvents)
        idx = x:size(b_bs_all{ne},2)+x-1;
        a=subplot(1,numPlots,ne);
        xaxis = linspace(histLen(ne,1),histLen(ne,2),size(b_bs_all{ne},2)); 
       
        imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(maxI,idx),[minC maxC]);
        colormap(flipud(red2blue));
        %colormap(bluewhitered);
        box off
        title(plotEvents{ne});
        if ne == numel(plotEvents)
            b=smallcolorbar(gca);
           b.Ticks = [minC 0 maxC];
            b.YLabel.String = 'Kernel value';
        end
        if ne ~= 1
            a.YTick = [];
        else
            ylabel('Neuron number (sorted by time of peak)')
        end
        x = x+size(b_bs_all{ne},2);
    end
end
hold on, plot(xaxis,ones(size(xaxis)).*find(maxI_orig(maxI)>size(temp2,2)/numel(plotEvents),1,'first'),'k') 
save(fullfile(fbasename,[savename, '.mat']), 'neuronIdx','pmat_all','b_all','cons','b_bs_all','recs');
%print(fullfile(whereAreWe('imaging'),'sex differences', [savename '_' norm '_' cat(2,cons{plotEvents}) '.pdf']),'-dpdf','-bestfit')
