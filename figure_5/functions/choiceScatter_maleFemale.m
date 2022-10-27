function choiceHist_maleFemale(recs_f,recs_m,frameRate,rawFlag,zscoreFlag,cohort)
if nargin < 6
    cohort = 'ACC_DMS_imaging';
end

a = .01;
bw = .05;
fbasename = fullfile(whereAreWe('imaging'),cohort);

savehere = fullfile(whereAreWe('figurecode'),'processed_data','fig5'); %generate save location
if ~isfolder(savehere)
    mkdir(savehere)
end

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));

histLength = [2 6];
[rewHists_m,nrewHists_m,rew_avg_m,nrew_avg_m,pvals_all_m] = choiceActivity(recs_m,frameRate,rawFlag,zscoreFlag,histLength,cohort);
[rewHists_f,nrewHists_f,rew_avg_f,nrew_avg_f,pvals_all_f] = choiceActivity(recs_f,frameRate,rawFlag,zscoreFlag,histLength,cohort);


save(fullfile(savehere,['choiceHist_maleFemale.mat']));




fbasename = fullfile(whereAreWe('imaging'));


savehere = fullfile(whereAreWe('figurecode'),'processed_data'); %generate save location
if ~isfolder(savehere)
    mkdir(savehere)
end

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); 

% load data 
load(fullfile(savehere,'outcomeHist_maleFemale_full.mat'));

%% Plot histogram of difference between reward and no reward on average
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;

screensize = get( groot, 'Screensize' );
screensize(4) = screensize(4)/2.5;
screensize(3) = screensize(3)/2.2;
screensize(2) = screensize(2)+screensize(4)-100;
figure('Position', screensize)
subplot(1,2,1), hold on

diff_avg_male = (rew_avg_m-nrew_avg_m)./(rew_avg_m+nrew_avg_m);
diff_avg_female = (rew_avg_f - nrew_avg_f)./(rew_avg_f+nrew_avg_f);

minD = min(cat(1,diff_avg_male,diff_avg_female));
maxD = max(cat(1,diff_avg_male,diff_avg_female));
if minD>maxD
    maxD = -minD;
elseif minD<maxD
    minD = -maxD;
end
minD = -1;
maxD = 1; 
hold on


[a,b] = histcounts(diff_avg_male,linspace(minD,maxD,20),'Normalization','probability');
[a2,b2] = histcounts(diff_avg_female,linspace(minD,maxD,20),'Normalization','probability');
plot(b(1:end-1),a,'Color',maleC);
plot(b(1:end-1),a2,'Color',femaleC);
plot(median(diff_avg_male),max(cat(2,a,a2))+.02,'v','Color',maleC,'MarkerSize',6)
plot(median(diff_avg_female),max(cat(2,a,a2))+.02,'v','Color',femaleC,'MarkerSize',6)

p_all = ranksum(diff_avg_male,diff_avg_female);
if p_all<.01
    title('All active neurons *')
else
    title('All active neurons')
end

xlabel('(Ipsi - Contra activity)/(Ipsi + Contra activity)')
ylabel('Proportion of neurons')

set(gca,'FontSize',16); 


end

function [rewHists,contraHists,ipsi_avg,contra_avg,pvals_all] = choiceActivity(recs,frameRate,rawFlag,zscoreFlag,histLength,cohort)

binSize = 1000/frameRate; % bin size for histogram
fbasename = fullfile(whereAreWe('imaging'),cohort);

% find significant outcome neurons
a = .01; % significance level
sigVer = sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_fig3_choice.mat',frameRate,1,zscoreFlag);
whichEvents = [4];


%% Time around event to average

histLengthF = histLength.*frameRate;

rewHists = [];
contraHists = [];
ipsi_avg = [];
contra_avg = [];
pvals_all = [];


for na = 1:numel(recs)
    load(fullfile(fbasename,recs{na}, ['behavEventsRegression_bs' num2str(binSize) '.mat']));
    info = matfile(fullfile(fbasename,recs{na},'info.mat'));
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
    
    %     % normalize gcamp
    %     if zscoreFlag
    %         mu = nanmean(dff,1);
    %         sigma = nanstd(dff,[],1);
    %         dff = (dff-repmat(mu,size(dff,1),1))./repmat(sigma, size(dff,1),1);
    %     end
    
    dff = dff(:,activeIdx);
    
    dff(dff<0) = 0;
    % ipsi lever press activity
    clear thishist
    nROI = size(dff,2);
    if info.ipsiSide == 1
        frames = tdtEvents.LLeverPress;
    elseif info.ipsiSide == 0
        frames = tdtEvents.RLeverPress;
    else
        keyboard
    end
    nevents = numel(frames);
    
    thishist = [];
    thisavg = [];
    for r = 1:nROI
        for t = 1:nevents
            thisidx = frames(t)-histLengthF(1):frames(t)+histLengthF(2);
            %thisidx = trialEnd(t)-frameRate*2:trialEnd(t);
            if thisidx(end) > length(dff) % sometimes last trial gets cut shorter than desired hist length
                thistrial = [dff(frames(t)-histLengthF(1):end,r) zeros(thisidx(end)-length(dff),1)];
            else thistrial = dff(thisidx,r);
            end
            thishist{t}(:,r) = thistrial;
            
            thisavg(r,t) = nanmean(thistrial);
        end
        
    end
    
    ipsiHists{na} = thishist;
    ipsi_avg = cat(1,ipsi_avg, nanmean(thisavg,2));
    ipsi = thisavg';
    
    % unrewarded outcome activity
    if info.ipsiSide == 1
        frames = tdtEvents.RLeverPress;
    elseif info.ipsiSide == 0
        frames = tdtEvents.LLeverPress;
    else
        keyboard
    end
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
    
    contraHists{na} = thishist;
    contra_avg = cat(1,contra_avg, nanmean(thisavg,2));
    contra = thisavg';
    
    clear temp
    %[~,temp] = ttest2(rew,nrew);
    for nn = 1:nROI
        [~,temp(nn)] = ttest2(ipsi(:,nn),contra(:,nn));
    end
    pvals_all = cat(1,pvals_all,temp');
    clear contra ipsi
    
end
end