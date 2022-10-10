function plotIndividualNeurons_fig2(recs,ver,frameRate,whichEvents,expt)

fbasename = fullfile(whereAreWe('imaging'),expt); 
binSize = 1000/frameRate;
fbasename_bs = fullfile(whereAreWe('bucket'),'DMS_Bandit', 'basis sets');
sigLevel = .05; 
%% Define events

rawFlag = 1;
for na = 1:numel(recs)
    savehere = fullfile(fbasename, 'individual ROIs', recs{na});
    idx = strfind(recs{na}, '\');
    if isempty(idx)
        idx = strfind(recs{na}, '/');
    end
    mkdir(fullfile(fbasename, 'individual ROIs', recs{na}(1:idx-1)));
    load(fullfile(fbasename,recs{na}, ['behavEventsRegression_bs' num2str(binSize) '.mat']));
    try dff = matfile(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'));
        dff = eval(sprintf('dff.dff_%d_raw1',frameRate));
    catch
        load(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'),'neuron');
        load(fullfile(fbasename,recs{na},'info.mat'),'droppedIdx','nFramesLog');
        info = matfile(fullfile(fbasename,recs{na},'info.mat'));
        frameRateO = info.frameRate;
        nframesLog = nFramesLog;
        if rawFlag %if rawFlag use the raw fluorescence
            dff = neuron.C_raw;
        else %otherwise use denoised fluorescence
            dff = neuron.C_raw;
        end
        
        if size(dff,2) > size(dff,1)
            dff = dff';
        end
        if droppedIdx == 0
            droppedIdx =[];
        end
        %if droppped frames had been removed for cnmfe, add back in NaN
        if numel(droppedIdx)>0 & numel(droppedIdx)+nframesLog > length(dff)
            temp = zeros(numel(droppedIdx)+nframesLog,size(dff,2));
            temp(droppedIdx,:) = NaN;
            temp(setdiff(1:numel(droppedIdx)+nframesLog,droppedIdx),:) = dff;
            dff = temp;
        elseif numel(droppedIdx) == 0 && nframesLog ~= length(dff)
            keyboard %this is probably T55 7/24 where I had to cut two frames to get inscopix to export tiff
        end
        %downsample if necessary
        if frameRate~=frameRateO
            dff = downsampleImaging(frameRateO, frameRate,dff);
        end
        eval(sprintf('dff_%d_raw1 = dff;', frameRate));
        save(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'), ['dff_' num2str(frameRate) '_raw1'], '-append');
    end
    
    % normalize gcamp
    
    mu = nanmean(dff,1);
    sigma = nanstd(dff,[],1);
    dff = (dff-repmat(mu,size(dff,1),1))./repmat(sigma, size(dff,1),1);
    
    
    % load parameters
    info = matfile(fullfile(fbasename,  recs{na}, 'info.mat'));
    frameRateO = info.frameRate;
    droppedIdx = info.droppedIdx;
    nframesLog = info.nframesLog;
    
    
    % only include neurons significantly modulated by at least one event (determined w/ F-stat on full v reduced regression)
    
    
    [cons, con_shift, time_back_orig, time_forward_orig, type1, shift_con, bsIDs,contVar,fullVar] = getEvents(ver,frameRate);
    
    % define events 
%Trial start
TrialStart = tdtEvents.trialStart;
%Nose poke
NPTimes = tdtEvents.nosePokeEntry;
NPTimesI = tdtEvents.nosePokeEntry(tdtEvents.choice==0);
NPTimesC = tdtEvents.nosePokeEntry(tdtEvents.choice==1);
NPTimesP = tdtEvents.nosePokeEntry(tdtEvents.reward==1);
NPTimesN = tdtEvents.nosePokeEntry(tdtEvents.reward==0&tdtEvents.choice~=-1);
idx = find(tdtEvents.reward==1)+1;
idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
NPPrevRew = tdtEvents.nosePokeEntry(idx);
idx = find(tdtEvents.reward==0&tdtEvents.choice~=-1)+1;
idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
NPPrevNRew = tdtEvents.nosePokeEntry(idx); 
%Lever presentation
LeverPresent = tdtEvents.leverPresentation;
LeverPresentI = tdtEvents.leverPresentation(tdtEvents.choice==0);
LeverPresentC = tdtEvents.leverPresentation(tdtEvents.choice==1);
LeverPresentP = tdtEvents.leverPresentation(tdtEvents.reward==1);
LeverPresentN = tdtEvents.leverPresentation(tdtEvents.reward==0&tdtEvents.choice~=-1);
%Lever press
LeverTimes = tdtEvents.leverPress;
LeverTimesI = tdtEvents.leverPress(tdtEvents.choice==0);
LeverTimesC = tdtEvents.leverPress(tdtEvents.choice==1);
LeverTimesP = tdtEvents.leverPress(tdtEvents.reward==1);
LeverTimesN = tdtEvents.leverPress(tdtEvents.reward==0&tdtEvents.choice~=-1);
%Outcome
CS = tdtEvents.outcome;
CSI = tdtEvents.outcome(tdtEvents.choice==0);
CSC = tdtEvents.outcome(tdtEvents.choice==1);
CS(tdtEvents.choice==-1) = []; %remove omitted trials
CSRew = tdtEvents.CSPlus;
CSNoRew = tdtEvents.CSMinus;

%Reward entry
RewardEnter = tdtEvents.rewardEntry;

%Reward exit
RewardExit = tdtEvents.rewardExit;
rewIdx = find(tdtEvents.reward==1);
rewIdx = setdiff(rewIdx,tdtEvents.trialExcludeReward); %in case of dropped frames/excluded trials
RewardExitI = tdtEvents.rewardExit(ismember(rewIdx,find(tdtEvents.choice==0)));
RewardExitC = tdtEvents.rewardExit(ismember(rewIdx,find(tdtEvents.choice==1)));


% load regression result
load(fullfile(fbasename,recs{na},['linReg_fs' num2str(frameRate) '_raw1_zscore1_basisSet_fig4_' ver '.mat']), 'b','con_iden','pvals')
load(fullfile(fbasename,recs{na}, ['M_linReg_fs' num2str(frameRate) '_raw1_zscore1_basisSet_fig4_' ver '.mat']), 'x_basic','event_times_mat');

pvals = cell2mat(pvals');
pmat = pvals<sigLevel;

for n = 1:size(event_times_mat,1)
    event_times_mat(n,event_times_mat(n,:)==1) = n;
end
    
    % generate response kernel
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
    
    
     % Generate kernel for each roi
     for nr = 1:size(dff,2)
         responseKernel{nr} = zeros(1,length([-time_back_orig*frameRate:(time_forward_orig+2)*frameRate]));
         for ne = 1:numel(cons)
             thisWeights = b{nr}(con_iden==ne);
             if iscell(bs)
                 tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*zscore(bs{ne}),2)';
             else
                 tempWeights = sum(repmat(thisWeights',size(bs,1),1).*zscore(bs),2)';
             end
             if con_shift(ne)
                 responseKernel{nr}((2*frameRate)+1:end) = responseKernel{nr}((2*frameRate)+1:end)+tempWeights;
             else
                 responseKernel{nr}(1:end-2*frameRate) = responseKernel{nr}(1:end-2*frameRate)+tempWeights;
             end
             individualKernel{nr}{ne} = tempWeights; 
         end
     end
    
    % Plot significant neurons for each event
    for ne = whichEvents
        appendFlag = 0;
        thisROI = find(pmat(:,ne)==1);
        if con_shift(ne)~=1
            histLengthF = [-time_back_orig time_forward_orig].*frameRate;
            %histLengthF = [-2 3].*frameRate;
        else
            histLengthF = [-time_back_orig+2 time_forward_orig+2].*frameRate;
            %histLengthF = [0 5].*frameRate;
        end
        
        if contains(cons{ne},'CS') || contains(cons{ne},'LeverTimesI') || contains(cons{ne},'LeverTimesC')
            if contains(cons{ne},'CS')
                frames1 = CSRew;
                frames2 = CSNoRew;
            else
                frames1 = LeverTimesI;
                frames2 = LeverTimesC;
            end
            
            for nr = thisROI'
                y_est = b{nr}'*cat(2,zeros(size(x_basic,1),1),zscore(x_basic))';
                thishist1 = [];
                thishist_est1 = [];
                thishist_event1 = [];
                thishist1_long = [];
                for t = 1:numel(frames1)
                    thisidx2 = frames1(t)+histLengthF(1)*2:frames1(t)+histLengthF(2)*2;
                    thisidx = frames1(t)+histLengthF(1):frames1(t)+histLengthF(2);
                    if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                        thistrial = [dff(frames1(t)-histLengthF(1):end,nr) zeros(thisidx(end)-length(dff),1)];
                    else thistrial = dff(thisidx,nr);
                    end
                    if thisidx2(end) > length(dff)
                        thistrial2 = [y_est(frames1(t)+histLengthF(1)*2:end) zeros(thisidx2(end)-length(dff),1)];
                        thistrial1_long = [dff(frames1(t)+histLengthF(1)*2:end,nr) zeros(thisidx2(end)-length(dff),1)];
                        
                        thistrial3 = [sum(event_times_mat(:,frames1(t)+histLengthF(1)*2:end)) zeros(thisidx2(end)-length(dff),1)];
                    else  thistrial2 = y_est(thisidx2); thistrial3 = sum(event_times_mat(:,thisidx2),1); thistrial1_long = dff(thisidx2,nr);
                    end
                    thishist1(t,:) = thistrial;
                    thishist1_long(t,:) = thistrial1_long;
                    thishist_est1(t,:) = thistrial2;
                    thishist_event1(t,:) = thistrial3;
                end
                
                thishist2 = [];
                thishist_est2 = [];
                thishist_event2 = [];
                thishist2_long = [];
                for t = 1:numel(frames2)
                    thisidx2 = frames2(t)+histLengthF(1)*2:frames2(t)+histLengthF(2)*2;
                    thisidx = frames2(t)+histLengthF(1):frames2(t)+histLengthF(2);
                    if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                        thistrial = [dff(frames2(t)-histLengthF(1):end,nr) zeros(thisidx(end)-length(dff),1)];
                    else thistrial = dff(thisidx,nr);
                    end
                    if thisidx2(end) > length(dff)
                        thistrial2 = [y_est(frames2(t)+histLengthF(1)*2:end) zeros(thisidx2(end)-length(dff),1)];
                        thistrial1_long = [dff(frames2(t)+histLengthF(1)*2:end,nr) zeros(thisidx2(end)-length(dff),1)];
                        
                        thistrial3 = [sum(event_times_mat(:,frames2(t)+histLengthF(1)*2:end)) zeros(thisidx2(end)-length(dff),1)];
                    else  thistrial2 = y_est(thisidx2); thistrial3 = sum(event_times_mat(:,thisidx2),1); thistrial1_long = dff(thisidx2,nr);
                    end
                    thishist2(t,:) = thistrial;
                    thishist2_long(t,:) = thistrial1_long;
                    thishist_est2(t,:) = thistrial2;
                    thishist_event2(t,:) = thistrial3;
                end
                
                
                
                nanidx = sum(isnan(thishist1),2)>1;
                thishist1(nanidx,:) = [];
                
                nanidx = sum(isnan(thishist2),2)>1;
                thishist2(nanidx,:) = [];
                
                maxC = max([max(max(thishist1)) abs(min(min(thishist1))) max(max(thishist2)) abs(min(min(thishist2)))]);
                minC = -maxC;
                xaxis = linspace(histLengthF(1)/frameRate,histLengthF(2)/frameRate,size(thishist1,2));
                xaxis2 = linspace(histLengthF(1)/frameRate,(histLengthF(2)/frameRate + 2),length(responseKernel{nr}));
                figure('Units','inches','Position',[5,5,6,2.5]),
                
                subplot(2,4,1), imagesc(xaxis,1:size(thishist1,1),thishist1, [-5 5]), colormap(flipud(red2blue)); smallcolorbar;
                title(cons{ne})
                subplot(2,4,5), imagesc(xaxis,1:size(thishist2,1),thishist2, [-5 5]), colormap(flipud(red2blue)); smallcolorbar;
                
                subplot(1,4,2), shadedErrorBar(xaxis,nanmean(thishist1),nanstd(thishist1)./sqrt(size(thishist1,1))),box off
                hold on, subplot(1,4,2), shadedErrorBar(xaxis,nanmean(thishist2),nanstd(thishist2)./sqrt(size(thishist2,1))),box off
                
                %subplot(1,4,3), plot(xaxis2, responseKernel{nr}), box off
                subplot(1,4,4), plot(xaxis,individualKernel{nr}{ne}), box off
                hold on, plot(xaxis,nanmean(thishist2))
                
                
                
                keyboard
                answer = inputdlg('check trials?');
                if answer{1} == '1'
                    idx = randi(size(thishist1,1),15,1);
                    figure()
                    thishist_event1(thishist_event1==0) = NaN;
                    for t =idx'
                        figure()
                        title(['trial' num2str(t)])
                        plot(thishist1_long(t,:),'k')
                        hold on
                        plot(thishist_est1(t,:),'b');
                        hold on
                        plot(thishist_event1(t,:),'ok');
                        pause()
                        close all
                    end
                    
                    idx = randi(size(thishist2,1),15,1);
                    figure()
                    thishist_event2(thishist_event2==0) = NaN;
                    for t =idx'
                        figure()
                        title(['trial' num2str(t)])
                        plot(thishist2_long(t,:),'k')
                        hold on
                        plot(thishist_est2(t,:),'b');
                        hold on
                        plot(thishist_event2(t,:),'ok');
                        pause()
                        close all
                    end
                end
                close all
            end
        else
            
            frames = eval(sprintf('%s',cons{ne}));
            nevents = numel(frames);
            for nr = thisROI'
                y_est = b{nr}'*cat(2,zeros(size(x_basic,1),1),zscore(x_basic))';
                
                
                thishist = [];
                thishist_est = [];
                thishist_event = [];
                thishist_long = [];
                for t = 1:nevents
                    
                        thisidx2 = frames(t)+histLengthF(1)*2:frames(t)+histLengthF(2)*2;
                        thisidx = frames(t)+histLengthF(1):frames(t)+histLengthF(2);
                        if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                            thistrial = [dff(frames(t)-histLengthF(1):end,nr) zeros(thisidx(end)-length(dff),1)];
                        else thistrial = dff(thisidx,nr);
                        end
                        if thisidx2(end) > length(dff)
                            thistrial2 = [y_est(frames(t)+histLengthF(1)*2:end) zeros(thisidx2(end)-length(dff),1)];
                            thistrial1_long = [dff(frames(t)+histLengthF(1)*2:end,nr) zeros(thisidx2(end)-length(dff),1)];
                            
                            thistrial3 = [sum(event_times_mat(:,frames(t)+histLengthF(1)*2:end)) zeros(thisidx2(end)-length(dff),1)];
                        else  thistrial2 = y_est(thisidx2); thistrial3 = sum(event_times_mat(:,thisidx2),1); thistrial1_long = dff(thisidx2,nr);
                        end
                        thishist(t,:) = thistrial;
                        thishist_long(t,:) = thistrial1_long;
                        thishist_est(t,:) = thistrial2;
                        thishist_event(t,:) = thistrial3;
                        
                    end
                    nanidx = sum(isnan(thishist),2)>1;
                    thishist(nanidx,:) = [];
                    maxC = max([max(max(thishist)) abs(min(min(thishist)))]);
                    minC = -maxC;
                    xaxis = linspace(histLengthF(1)/frameRate,histLengthF(2)/frameRate,size(thishist,2));
                    xaxis2 = linspace(histLengthF(1)/frameRate,(histLengthF(2)/frameRate + 2),length(responseKernel{nr}));
                    figure('Units','inches','Position',[5,5,6,2.5]), subplot(1,4,1), imagesc(xaxis,1:size(thishist,1),thishist, [-5 5]), colormap(flipud(red2blue)); smallcolorbar;
                    subplot(1,4,2), shadedErrorBar(xaxis,nanmean(thishist),nanstd(thishist)./sqrt(size(thishist,1))),box off
                    %subplot(1,4,3), plot(xaxis2, responseKernel{nr}), box off
                    subplot(1,4,4), plot(xaxis,individualKernel{nr}{ne}), box off
                    hold on, plot(xaxis,nanmean(thishist))
                    title(cons{ne})
                    %            if appendFlag==0
                    %                export_fig([savehere, cons{ne}, '.pdf'])
                    %                appendFlag = 1;
                    %            else
                    %                export_fig([savehere, cons{ne}, '.pdf'],'-append')
                    %            end
                    keyboard
                    answer = inputdlg('check trials?');
                    if answer{1} == '1'
                        idx = randi(size(thishist,1),15,1);
                        figure()
                        thishist_event(thishist_event==0) = NaN;
                        for t =idx'
                            figure()
                            title(['trial' num2str(t)])
                            plot(thishist_long(t,:),'k')
                            hold on
                            plot(thishist_est(t,:),'b');
                            hold on
                            plot(thishist_event(t,:),'ok');
                            pause()
                            close all
                        end
                    end
                    close all
                end
            end
        end
        
    
    
%     % plot activity divided by choice
%     events = {'trialStart';'nosePokeEntry';'leverPresentation';'leverPress';'outcome'};
%     histLength = [3 2; 2 2; 2 2; 2 2; 0 3];
%     
%     histLengthF = histLength.*frameRate;
%     for r = 1:nROI
%         sigMat = [0, pmat(r,1), pmat(r,2), pmat(r,3)| pmat(r,4), pmat(r,5)|pmat(r,6)|pmat(r,7)];
%         screensize = get( groot, 'Screensize' );
%         figure('Position', screensize)
%         for ne = 1:numel(events)
%             thishist_L = [];
%             thishist_R = [];
%             thisavg_L = [];
%             thisavg_R = [];
%             
%             frames = eval(sprintf('tdtEvents.%s(tdtEvents.choice==1)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_L(t,:) = thistrial;
%             end
%             frames = eval(sprintf('tdtEvents.%s(tdtEvents.choice==0)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_R(t,:) = thistrial;
%             end
%             clear p p2
%             subplot(3, numel(events),ne)
%             p=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_R,2)), nanmean(thishist_R,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_R,2)), nanmean(thishist_R,1)+nanstd(thishist_R,1)./sqrt(size(thishist_R,1)), ':', 'Color', p.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_R,2)), nanmean(thishist_R,1)-nanstd(thishist_R,1)./sqrt(size(thishist_R,1)), ':', 'Color', p.Color)
%             p2=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_L,2)), nanmean(thishist_L,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_L,2)),nanmean(thishist_L,1)+nanstd(thishist_L,1)./sqrt(size(thishist_L,1)), ':', 'Color', p2.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_L,2)), nanmean(thishist_L,1)-nanstd(thishist_L,1)./sqrt(size(thishist_L,1)), ':', 'Color', p2.Color)
%             plot([0 0], [-1, 1],'k')
%             if ne == numel(events)
%                 legend([p p2], {'Right';'Left'})
%             end
%             if ne == 1
%                 ylabel('zscore gcamp')
%             end
%             
%             if sigMat(ne) == 1
%                 title([events{ne} '*']);
%             else
%                 title(events{ne})
%             end
%             
%             subplot(3,numel(events),ne+numel(events))
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_R,2)),1:size(thishist_R,1),thishist_R,[min(min([thishist_R; thishist_L])) max(max([thishist_R; thishist_L]))]);
%             
%             if ne == 1
%                 ylabel(sprintf('Right trials \n trial number'))
%             end
%             subplot(3,numel(events),ne+numel(events)*2)
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_L,2)),1:size(thishist_L,1),thishist_L,[min(min([thishist_R; thishist_L])) max(max([thishist_R; thishist_L]))]);
%             
%             if ne == 1
%                 ylabel(sprintf('Left trials \n trial number'))
%             end
%             smallcolorbar(gca,'eastoutside');
%             xlabel('Time (s)')
%             
%         end
%         if r == 1
%             export_fig(savehere,'-pdf')
%         else
%             export_fig(savehere,'-pdf','-append')        
%         end
%         
%         
%         
%         screensize = get( groot, 'Screensize' );
%         figure('Position', screensize)
%         for ne = 1:numel(events)
%            
%              thishist_rew = [];
%             thishist_nrew = [];
%             
%             frames = eval(sprintf('tdtEvents.%s(tdtEvents.reward==1)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_rew(t,:) = thistrial;
%             end
%             frames = eval(sprintf('tdtEvents.%s(tdtEvents.reward==0&tdtEvents.choice~=-1)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_nrew(t,:) = thistrial;
%             end
%             subplot(3, numel(events),ne)
%             clear p p2
%             p=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)), nanmean(thishist_rew,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)), nanmean(thishist_rew,1)+nanstd(thishist_rew,1)./sqrt(size(thishist_rew,1)), ':', 'Color', p.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)), nanmean(thishist_rew,1)-nanstd(thishist_rew,1)./sqrt(size(thishist_rew,1)), ':', 'Color', p.Color)
%             p2=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)), nanmean(thishist_nrew,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)),nanmean(thishist_nrew,1)+nanstd(thishist_nrew,1)./sqrt(size(thishist_nrew,1)), ':', 'Color', p2.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)), nanmean(thishist_nrew,1)-nanstd(thishist_nrew,1)./sqrt(size(thishist_nrew,1)), ':', 'Color', p2.Color)
%             plot([0 0], [-1, 1],'k')
%             if ne == numel(events)
%                 legend([p p2], {'Reward';'No reward'})
%             end
%             if ne == 1
%                 ylabel('zscore gcamp')
%             end
%             
%             if sigMat(ne) == 1
%                 title([events{ne} '*']);
%             else
%                 title(events{ne})
%             end
%             
%             subplot(3,numel(events),ne+numel(events))
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)),1:size(thishist_rew,1),thishist_rew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
%             
%             if ne == 1
%                 ylabel(sprintf('Rewarded trials \n trial number'))
%             end
%             subplot(3,numel(events),ne+numel(events)*2)
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)),1:size(thishist_nrew,1),thishist_nrew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
%             
%             if ne == 1
%                 ylabel(sprintf('Unrewarded trials \n trial number'))
%             end
%             smallcolorbar(gca,'eastoutside');
%             xlabel('Time (s)')
%         end
%         export_fig(savehere,'-pdf','-append')
%         
%         
%         % plot activity divided by outcome and choice
%         clear thishit_L thishistR
%         screensize = get( groot, 'Screensize' );
%         figure('Position', screensize)
%         for ne = 1:numel(events)
%             thishist_rewL = [];
%             thishist_rewR = [];
%             thishist_nrewL = [];
%             thishist_nrewR = [];
%             idx = find(tdtEvents.reward==1);
%             idx = idx(ismember(idx,find(tdtEvents.choice==1))); %left and rewarded
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_rewL(t,:) = thistrial;
%             end
%             idx = find(tdtEvents.reward==1);
%             idx = idx(ismember(idx,find(tdtEvents.choice==0))); %right and rewarded
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_rewR(t,:) = thistrial;
%             end
%             %left and unrewarded
%             idx = find(tdtEvents.reward==0);
%             idx = idx(ismember(idx,find(tdtEvents.choice==1))); 
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_nrewL(t,:) = thistrial;
%             end
%             %left and unrewarded
%             idx = find(tdtEvents.reward==0);
%             idx = idx(ismember(idx,find(tdtEvents.choice==0))); 
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_nrewR(t,:) = thistrial;
%             end
%             
%             subplot(3, numel(events),ne)
%             clear p p2 p3 p4
%             p=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1)+nanstd(thishist_rewL,1)./sqrt(size(thishist_rewL,1)), ':', 'Color', p.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1)-nanstd(thishist_rewL,1)./sqrt(size(thishist_rewL,1)), ':', 'Color', p.Color)
%             p2=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)), nanmean(thishist_rewR,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)),nanmean(thishist_rewR,1)+nanstd(thishist_rewR,1)./sqrt(size(thishist_rewR,1)), ':', 'Color', p2.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)), nanmean(thishist_rewR,1)-nanstd(thishist_rewR,1)./sqrt(size(thishist_rewR,1)), ':', 'Color', p2.Color)
%             p3=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1)+nanstd(thishist_nrewL,1)./sqrt(size(thishist_nrewL,1)), ':', 'Color', p3.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1)-nanstd(thishist_nrewL,1)./sqrt(size(thishist_nrewL,1)), ':', 'Color', p3.Color)
%             p4=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)), nanmean(thishist_nrewR,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)),nanmean(thishist_nrewR,1)+nanstd(thishist_nrewR,1)./sqrt(size(thishist_nrewR,1)), ':', 'Color', p4.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)), nanmean(thishist_nrewR,1)-nanstd(thishist_nrewR,1)./sqrt(size(thishist_nrewR,1)), ':', 'Color', p4.Color)
%           
%             
%             plot([0 0], [-1, 1],'k')
%             if ne == numel(events)
%                 l=legend([p p2 p3 p4], {'Rewarded left choice';'Rewarded right choice'; 'Unrewarded left choice';'Unrewarded right choice'});
%                 l.Box = 'off';
%             end
%             if ne == 1
%                 ylabel('zscore gcamp')
%             end
%             
%             if sigMat(ne) == 1
%                 title([events{ne} '*']);
%             else
%                 title(events{ne})
%             end
%             
% %             subplot(3,numel(events),ne+numel(events))
% %             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)),1:size(thishist_rew,1),thishist_rew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
% %             
% %             if ne == 1
% %                 ylabel(sprintf('Rewarded trials \n trial number'))
% %             end
% %             subplot(3,numel(events),ne+numel(events)*2)
% %             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)),1:size(thishist_nrew,1),thishist_nrew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
% %             
% %             if ne == 1
% %                 ylabel(sprintf('Unrewarded trials \n trial number'))
% %             end
% %             smallcolorbar(gca,'eastoutside');
% %             xlabel('Time (s)')
%         end
%         export_fig(savehere,'-pdf','-append')
%         
%         
%         % plot activity divided by previos outcome
%         
%         
%         screensize = get( groot, 'Screensize' );
%         figure('Position', screensize)
%         for ne = 1:numel(events)
%             thishist_L = [];
%             thishist_R = [];
%             thishist_prevRew = [];
%             thishist_prevNrew = [];
%             
%             idx = find(tdtEvents.reward==1)+1;
%             idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             frames(isnan(frames)) = [];
%             nevents = numel(frames);
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_prevRew(t,:) = thistrial;
%             end
%             idx = find(tdtEvents.reward==0&tdtEvents.choice~=-1)+1;
%             idx = idx(ismember(idx,1:numel(tdtEvents.reward)));
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             frames(isnan(frames)) = [];
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_prevNrew(t,:) = thistrial;
%             end
%             subplot(3, numel(events),ne)
%             clear p p2
%              p=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevRew,2)), nanmean(thishist_prevRew,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevRew,2)),nanmean(thishist_prevRew,1)+nanstd(thishist_prevRew,1)./sqrt(size(thishist_prevRew,1)), ':', 'Color', p.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevRew,2)), nanmean(thishist_prevRew,1)-nanstd(thishist_prevRew,1)./sqrt(size(thishist_prevRew,1)), ':', 'Color', p.Color)
%              p2=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevNrew,2)), nanmean(thishist_prevNrew,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevNrew,2)), nanmean(thishist_prevNrew,1)+nanstd(thishist_prevNrew,1)./sqrt(size(thishist_prevNrew,1)), ':', 'Color', p2.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevNrew,2)), nanmean(thishist_prevNrew,1)-nanstd(thishist_prevNrew,1)./sqrt(size(thishist_prevNrew,1)), ':', 'Color', p2.Color)
%           
%             plot([0 0], [-1, 1],'k')
%             if ne == 1
%                 ylabel('zscore gcamp')
%             end
%             if ne == numel(events)
%                 legend([p p2], {'Previous reward';'Previous no reward'})
%             end
%             if sigMat(ne) == 1
%                 title([events{ne} '*']);
%             else
%                 title(events{ne})
%             end
%             
%             subplot(3,numel(events),ne+numel(events)*2)
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevNrew,2)),1:size(thishist_prevNrew,1),thishist_prevNrew,[min(min([thishist_prevNrew; thishist_prevRew])) max(max([thishist_prevNrew; thishist_prevRew]))]);
%             
%             if ne == 1
%                 ylabel(sprintf('Previously nrewarded trials \n trial number'))
%             end
%             subplot(3,numel(events),ne+numel(events))
%             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_prevRew,2)),1:size(thishist_prevRew,1),thishist_prevRew,[min(min([thishist_prevNrew; thishist_prevRew])) max(max([thishist_prevNrew; thishist_prevRew]))]);
%             if ne == 1
%                 ylabel(sprintf('Previously rewarded trials \n trial number'))
%             end
%             
%             smallcolorbar(gca,'eastoutside');
%             xlabel('Time (s)')
%         end
%         export_fig(savehere,'-pdf','-append')
%         
%         
%          % plot activity divided by previous outcome and choice
%         clear thishit_L thishistR
%         screensize = get( groot, 'Screensize' );
%         figure('Position', screensize)
%         for ne = 1:numel(events)
%             thishist_rewL = [];
%             thishist_rewR = [];
%             thishist_nrewL = [];
%             thishist_nrewR = [];
%             idx = find(tdtEvents.reward==1)+1;
%             idx = idx(ismember(idx,find(tdtEvents.choice==1))); %left and rewarded
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_rewL(t,:) = thistrial;
%             end
%             idx = find(tdtEvents.reward==1)+1;
%             idx = idx(ismember(idx,find(tdtEvents.choice==0))); %right and rewarded
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_rewR(t,:) = thistrial;
%             end
%             %left and unrewarded
%             idx = find(tdtEvents.reward==0)+1;
%             idx = idx(ismember(idx,find(tdtEvents.choice==1))); 
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_nrewL(t,:) = thistrial;
%             end
%             %left and unrewarded
%             idx = find(tdtEvents.reward==0)+1;
%             idx = idx(ismember(idx,find(tdtEvents.choice==0))); 
%             frames = eval(sprintf('tdtEvents.%s(idx)',events{ne}));
%             nevents = numel(frames);
%             
%             for t = 1:nevents
%                 thisidx = frames(t)-histLengthF(ne,1):frames(t)+histLengthF(ne,2);
%                 if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
%                     thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
%                 else thistrial = dff(thisidx,r);
%                 end
%                 thishist_nrewR(t,:) = thistrial;
%             end
%             
%             subplot(3, numel(events),ne)
%             clear p p2 p3 p4
%             p=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1)+nanstd(thishist_rewL,1)./sqrt(size(thishist_rewL,1)), ':', 'Color', p.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewL,2)), nanmean(thishist_rewL,1)-nanstd(thishist_rewL,1)./sqrt(size(thishist_rewL,1)), ':', 'Color', p.Color)
%             p2=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)), nanmean(thishist_rewR,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)),nanmean(thishist_rewR,1)+nanstd(thishist_rewR,1)./sqrt(size(thishist_rewR,1)), ':', 'Color', p2.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rewR,2)), nanmean(thishist_rewR,1)-nanstd(thishist_rewR,1)./sqrt(size(thishist_rewR,1)), ':', 'Color', p2.Color)
%             p3=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1)+nanstd(thishist_nrewL,1)./sqrt(size(thishist_nrewL,1)), ':', 'Color', p3.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewL,2)), nanmean(thishist_nrewL,1)-nanstd(thishist_nrewL,1)./sqrt(size(thishist_nrewL,1)), ':', 'Color', p3.Color)
%             p4=plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)), nanmean(thishist_nrewR,1));
%             hold on
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)),nanmean(thishist_nrewR,1)+nanstd(thishist_nrewR,1)./sqrt(size(thishist_nrewR,1)), ':', 'Color', p4.Color)
%             plot(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrewR,2)), nanmean(thishist_nrewR,1)-nanstd(thishist_nrewR,1)./sqrt(size(thishist_nrewR,1)), ':', 'Color', p4.Color)
%           
%             
%             plot([0 0], [-1, 1],'k')
%             if ne == numel(events)
%                 l=legend([p p2 p3 p4], {'Prev Rew & left choice';' Prev Rew & right choice'; 'Prev Unrew & left choice';'Prev Unrew & right choice'});
%                 l.Box = 'off';
%             end
%             if ne == 1
%                 ylabel('zscore gcamp')
%             end
%             
%             if sigMat(ne) == 1
%                 title([events{ne} '*']);
%             else
%                 title(events{ne})
%             end
%             
% %             subplot(3,numel(events),ne+numel(events))
% %             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_rew,2)),1:size(thishist_rew,1),thishist_rew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
% %             
% %             if ne == 1
% %                 ylabel(sprintf('Rewarded trials \n trial number'))
% %             end
% %             subplot(3,numel(events),ne+numel(events)*2)
% %             imagesc(linspace(-histLength(ne,1),histLength(ne,2),size(thishist_nrew,2)),1:size(thishist_nrew,1),thishist_nrew,[min(min([thishist_rew; thishist_nrew])) max(max([thishist_rew; thishist_nrew]))]);
% %             
% %             if ne == 1
% %                 ylabel(sprintf('Unrewarded trials \n trial number'))
% %             end
% %             smallcolorbar(gca,'eastoutside');
% %             xlabel('Time (s)')
%         end
%         export_fig(savehere,'-pdf','-append')
%         
%         
%         close all
%     end
%     
end     
        
      