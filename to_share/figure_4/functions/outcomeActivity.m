function [rewHists,nrewHists,rew_avg,nrew_avg,pvals] = outcomeActivity(recs,frameRate,rawFlag,histLength,fbasename)

binSize = 1000/frameRate; % bin size for histogram 


%% Time around event to average

histLengthF = histLength.*frameRate;

rewHists = [];
nrewHists = [];
rew_avg = [];
nrew_avg = [];
pvals = [];

for na = 1:numel(recs)
    load(fullfile(fbasename,recs{na}, ['behavEventsRegression_bs' num2str(binSize) '.mat']));

    % load parameters
    info = matfile(fullfile(fbasename,  recs{na}, 'info.mat'));
    frameRateO = info.frameRate;
    droppedIdx = info.droppedIdx;
    nframesLog = info.nframesLog;
    
    % load indices of active neurons
    load(fullfile(fbasename, recs{na}, 'activeNeurons.mat'))
    
    try dff = matfile(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'));
        dff = eval(sprintf('dff.dff_%s_raw%s',num2str(frameRate), num2str(rawFlag)));
    catch
        load(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'),'neuron');
        dff = neuron.C_raw;
        
        
        if size(dff,2) > size(dff,1)
            dff = dff';
        end
        info = matfile(fullfile(fbasename,recs{na},'info.mat'));
        droppedIdx = info.droppedIdx;
        nframesLog = info.nframesLog;
        frameRateO = info.frameRate;
        %if droppped frames had been removed for cnmfe, add back in NaN
        clear temp
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
            dff = downsampleImaging(frameRateO, frameRate,dff);
        end
        eval(sprintf('dff_%d_raw1 = dff;', frameRate));
        save(fullfile(fbasename,recs{na},'analyzed_neuron_final.mat'), ['dff_' num2str(frameRate) '_raw' num2str(rawFlag)], '-append');
    end

    
    dff = dff(:,activeIdx);
    % full trial
    if histLengthF == 0
        % rewarded outcome activity
        clear thishist
        nROI = size(dff,2);
        frames = tdtEvents.rewardDelivery;
        nevents = numel(frames);
        %find end of the ITI (defined as nose poke of next trial)
        try
            trialEnd = tdtEvents.nosePokeEntry(find(tdtEvents.reward==1)+1);
        catch
            trialEnd = tdtEvents.nosePokeEntry(find(tdtEvents.reward(1:end-1)==1)+1);
            trialEnd = cat(2,trialEnd,size(dff,1));
        end
        thishist = [];
        thisavg = [];
        for r = 1:nROI
            for t = 1:nevents
                thisidx = frames(t)-histLengthF(1):trialEnd(t);
                if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                    thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
                else thistrial = dff(thisidx,r);
                end
                thishist{t}(:,r) = thistrial;
               
                thisavg(r,t) = nanmean(thistrial);
            end
           
        end
        
        rewHists{na} = thishist;
        rew_avg = cat(1,rew_avg, nanmean(thisavg,2));
        rew = thisavg';
        % unrewarded outcome activity
        frames = tdtEvents.CSMinus;
        nevents = numel(frames);
        try
            trialEnd = tdtEvents.nosePokeEntry(find(tdtEvents.reward==0&tdtEvents.choice~=-1)+1);
        catch
            trialEnd = tdtEvents.nosePokeEntry(find(tdtEvents.reward(1:end-1)==0&tdtEvents.choice(1:end-1)~=-1)+1);
            trialEnd = cat(2,trialEnd,size(dff,1));
        end
        thishist = [];
        thisavg = [];
     
        for r = 1:nROI
            for t = 1:nevents
                thisidx = frames(t)-histLengthF(1):trialEnd(t);                
                if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                    thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
                else thistrial = dff(thisidx,r);
                end
                thishist{t}(:,r) = thistrial;
             
                thisavg(r,t) = nanmean(thistrial);
            end
            
        end
        
        nrewHists{na} = thishist;
        nrew_avg = cat(1,nrew_avg, nanmean(thisavg,2));
        nrew = thisavg';

        clear temp
        for nn = 1:nROI
            [~,temp(nn)] = ttest2(rew(:,nn),nrew(:,nn));
        end
        pvals = cat(1,pvals,temp');
        clear rew nrew
    else
        % rewarded outcome activity
        nROI = size(dff,2);
        frames = tdtEvents.rewardDelivery;
        nevents = numel(frames);
        
        thishist = [];
        thisavg = [];
       
        for r = 1:nROI
            for t = 1:nevents
                thisidx = frames(t)-histLengthF(1):frames(t)+histLengthF(2);
                if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                    thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
                else thistrial = dff(thisidx,r);
                end
                thishist{t}(:,r) = thistrial;
               
                thisavg(r,t) = nanmean(thistrial);
            end
          
        end
        
        rewHists{na} = thishist;
        rew_avg = cat(1,rew_avg, nanmean(thisavg,2));
        rew = thisavg';
        
        % unrewarded outcome activity
        frames = tdtEvents.CSMinus;
        nevents = numel(frames);
       
        thishist = [];
        thisavg = [];
        ct = 1;
        thishist_outcome = [];
        for r = 1:nROI
            for t = 1:nevents
                thisidx = frames(t)-histLengthF(1):frames(t)+histLengthF(2);                
                if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                    thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
                else thistrial = dff(thisidx,r);
                end
                thishist{t}(:,r) = thistrial;
                
                thisavg(r,t) = nanmean(thistrial);
            end
            
           
        end
        
        nrewHists{na} = thishist;
        nrew_avg = cat(1,nrew_avg, nanmean(thisavg,2));
        nrew = thisavg';
        
        clear temp
        for nn = 1:nROI
            [~,temp(nn)] = ttest2(rew(:,nn),nrew(:,nn));
        end
        pvals = cat(1,pvals,temp');
        clear rew nrew
    end
end
end