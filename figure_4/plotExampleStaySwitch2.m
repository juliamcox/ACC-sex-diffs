function plotExampleStaySwitch2(out,eventDur,plotC)
for ne = 1:size(eventDur,1)
    nIdx = find(~logical(out.negKernel{1}));
    nIdx = [1:numel(out.negKernel{1})]';
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out.dff_rew_willStay{ne,1},2));
    
    for n = nIdx'
         figure(); hold on
        shadedErrorBar(xaxis,nanmean(out.dff_nrew_willStay_trials{ne}{n}),nanstd(out.dff_nrew_willStay_trials{ne}{n})./sqrt(size(out.dff_nrew_willStay_trials{ne}{ne},1)),'LineProps',{'--'});
        shadedErrorBar(xaxis,nanmean(out.dff_nrew_willSwitch_trials{ne}{n}),nanstd(out.dff_nrew_willSwitch_trials{ne}{n})./sqrt(size(out.dff_nrew_willStay_trials{ne}{ne},1)));
     
        
        prevReward = [out.reward{n}(2:end) NaN];
        thisTrials_rew = out.thistrials_orig{n}(out.reward{n}==1,:);
        thisTrials_nrew = out.thistrials_orig{n}(out.reward{n}==0,:);
        thisTrials_switch = out.thistrials_orig{n}(out.willStay{n}==0,:);
        thisTrials_stay = out.thistrials_orig{n}(out.willStay{n}==1,:);
        thisTrials_switchRew = out.thistrials_orig{n}(out.willStay{n}==0&out.reward{n}==1,:);
        thisTrials_stayRew = out.thistrials_orig{n}(out.willStay{n}==1&out.reward{n}==1,:);
        
        thisTrials_switchNRew = out.thistrials_orig{n}(out.willStay{n}==0&out.reward{n}==0,:);
        thisTrials_stayNRew = out.thistrials_orig{n}(out.willStay{n}==1&out.reward{n}==0,:);
        
        trialStart = out.trialStart{n};
        trialStart = [trialStart(2:end) NaN];
        nosePoke = out.nosePokeEntry{n};
        nosePoke = [nosePoke(2:end) NaN];
        figure('Position',[440 46 560 752]);
        % Plot reward
        subplot(4,2,1); hold on
        thisTimestamp = out.latency{n}(out.reward{n}==1);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_rew,1),thisTrials_rew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.reward{n}==1)-out.outcomeTime{n}(out.reward{n}==1);
                thisTimestamp(thisTimestamp>8) = NaN;
                plot(thisTimestamp(sortIdx),1:size(thisTrials_rew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        thisTimestamp = thisTimestamp + out.latency{n}(out.reward{n}==1)';
%         thisTimestamp = -out.latency{n}(out.reward{n}==1)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
        thisTimestamp(thisTimestamp>8) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_rew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('rewarded trials')
        axis tight
        
        % Plot no reward
        subplot(4,2,2); hold on
        thisTimestamp = out.latency{n}(out.reward{n}==0);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_nrew,1),thisTrials_nrew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
        thisTimestamp = trialStart(out.reward{n}==0)-out.outcomeTime{n}(out.reward{n}==0);
        thisTimestamp(thisTimestamp>8) = NaN;
        %         plot(thisTimestamp(sortIdx),1:size(thisTrials_nrew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
        thisTimestamp = thisTimestamp + out.latency{n}(out.reward{n}==0)';
        thisTimestamp(thisTimestamp>8) = NaN;
        
        %         thisTimestamp = -out.latency{n}(out.reward{n}==0)';
        %         thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_nrew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('unrewarded trials')
        axis tight
        
        
        
        % stay trials
        subplot(4,2,3); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==1);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_stay,1),thisTrials_stay(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==1)-out.outcomeTime{n}(out.willStay{n}==1);
                thisTimestamp(thisTimestamp>8) = NaN;
          %      plot(thisTimestamp(sortIdx),1:size(thisTrials_stay,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
                thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==1)';
                thisTimestamp(thisTimestamp>8) = NaN;
%         thisTimestamp = -out.latency{n}(out.willStay{n}==1)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_stay,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('stay trials')
        
        subplot(4,2,4); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==0);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_switch,1),thisTrials_switch(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==0)-out.outcomeTime{n}(out.willStay{n}==0);
                thisTimestamp(thisTimestamp>8) = NaN;
           %     plot(thisTimestamp(sortIdx),1:size(thisTrials_switch,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
                thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==0)';
                thisTimestamp(thisTimestamp>8) = NaN;
%         thisTimestamp = -out.latency{n}(out.willStay{n}==0)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_switch,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('switch trials')
        
        % stay rewarded trials
        subplot(4,2,5); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==1&out.reward{n}==1);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_stayRew,1),thisTrials_stayRew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==1&out.reward{n}==1)-out.outcomeTime{n}(out.willStay{n}==1&out.reward{n}==1);
                thisTimestamp(thisTimestamp>8) = NaN;
        %         plot(thisTimestamp(sortIdx),1:size(thisTrials_stayRew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
                thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==1&out.reward{n}==1)';
                thisTimestamp(thisTimestamp>8) = NaN;
       % thisTimestamp = -out.latency{n}(out.willStay{n}==1&out.reward{n}==1)';
        %thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_stayRew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('rewarded stay trials')
        
        % switch rewarded trials
        subplot(4,2,6); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==0&out.reward{n}==1);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_switchRew,1),thisTrials_switchRew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==0&out.reward{n}==1)-out.outcomeTime{n}(out.willStay{n}==0&out.reward{n}==1);
                thisTimestamp(thisTimestamp>8) = NaN;
        %   plot(thisTimestamp(sortIdx),1:size(thisTrials_switchRew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
                thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==0&out.reward{n}==1)';
                thisTimestamp(thisTimestamp>8) = NaN;
%         thisTimestamp = -out.latency{n}(out.willStay{n}==0&out.reward{n}==1)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_switchRew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('rewarded switch trials')
        
        subplot(4,2,7); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==1&out.reward{n}==0);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_stayNRew,1),thisTrials_stayNRew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==1&out.reward{n}==0)-out.outcomeTime{n}(out.willStay{n}==1&out.reward{n}==0);
                thisTimestamp(thisTimestamp>8) = NaN;
       % plot(thisTimestamp(sortIdx),1:size(thisTrials_stayNRew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
                thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==1&out.reward{n}==0)';
                thisTimestamp(thisTimestamp>8) = NaN;
%         thisTimestamp = -out.latency{n}(out.willStay{n}==1&out.reward{n}==0)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
        plot(thisTimestamp(sortIdx),1:size(thisTrials_stayNRew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('unrewarded stay trials')
        
        % switch rewarded trials
        subplot(4,2,8); hold on
        thisTimestamp = out.latency{n}(out.willStay{n}==0&out.reward{n}==0);
        [~,sortIdx] = sort(thisTimestamp);
        sortIdx = 1:numel(thisTimestamp);
        imagesc(xaxis,1:size(thisTrials_switchNRew,1),thisTrials_switchNRew(sortIdx,:),[0 8])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
                thisTimestamp = trialStart(out.willStay{n}==0&out.reward{n}==0)-out.outcomeTime{n}(out.willStay{n}==0&out.reward{n}==0);
                thisTimestamp(thisTimestamp>8) = NaN;
        %         plot(thisTimestamp(sortIdx),1:size(thisTrials_switchNRew,1),'sm','MarkerFaceColor','m','MarkerSize',2)
        axis tight
        thisTimestamp = thisTimestamp + out.latency{n}(out.willStay{n}==0&out.reward{n}==0)';
        thisTimestamp(thisTimestamp>8) = NaN;
%         thisTimestamp = -out.latency{n}(out.willStay{n}==0&out.reward{n}==0)';
%         thisTimestamp(thisTimestamp<-4) = NaN;
    %    plot(thisTimestamp(sortIdx),1:size(thisTrials_switchNRew,1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('unrewarded switch trials')
        
        pause()
        close all
        
    end
    
    
    
    
    
end