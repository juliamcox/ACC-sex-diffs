function plotExampleStaySwitch(out,eventDur,plotC)
for ne = 1:size(eventDur,1)
    nIdx = find(~logical(out.negKernel{1}));
    nIdx = [1:numel(out.negKernel{1})]'; 
    xaxis = linspace(-eventDur(ne,1),eventDur(ne,2),size(out.dff_rew_willStay{ne,1},2));
    
    for n = nIdx'
        figure(); hold on
        shadedErrorBar(xaxis,nanmean(out.dff_nrew_willStay_trials{ne}{n}),nanstd(out.dff_nrew_willStay_trials{ne}{n})./sqrt(size(out.dff_nrew_willStay_trials{ne}{ne},1)),'LineProps',{'--'});
        shadedErrorBar(xaxis,nanmean(out.dff_nrew_willSwitch_trials{ne}{n}),nanstd(out.dff_nrew_willSwitch_trials{ne}{n})./sqrt(size(out.dff_nrew_willStay_trials{ne}{ne},1)));
        
        thisTrials = cat(1,out.dff_nrew_willSwitch_trials{ne}{n}, out.dff_nrew_willStay_trials{ne}{n});
        maxC = max(max(thisTrials));
        minC = min(min(thisTrials));
        if maxC >8
            maxC=8;
        end
%         figure()
%         subplot(2,1,1)
%         imagesc(xaxis,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),out.dff_nrew_willSwitch_trials{ne}{n},[minC,maxC])
%         colormap(flipud(red2blue))
%         smallcolorbar()
%         thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==0);
%         hold on;
%         %plot(thisTimestamp+3,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
%         title('switch trials')
%         subplot(2,1,2)
%         imagesc(xaxis,1:size(out.dff_nrew_willStay_trials{ne}{n},1),out.dff_nrew_willStay_trials{ne}{n},[minC,maxC])
%         colormap(flipud(red2blue))
%         title('stay trials')
%         hold on
%         thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==1);
%         hold on;
%        % plot(thisTimestamp+3,1:size(out.dff_nrew_willStay_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
%         smallcolorbar()
%         
         
        figure();
        subplot(2,1,1)
        thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==0);
        thistrial = out.thistrials_resid{n}(out.reward{n}==0&out.willStay{n}==0,:);
        [~,sortIdx] = sort(thisTimestamp);
        imagesc(xaxis,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),thistrial(sortIdx,:),[minC,maxC])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
        plot(thisTimestamp(sortIdx)+3,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('switch trials')
        subplot(2,1,2)
        thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==1);
        thistrial = out.thistrials_resid{n}(out.reward{n}==0&out.willStay{n}==1,:);
        
        [~,sortIdx] = sort(thisTimestamp);
        imagesc(xaxis,1:size(out.dff_nrew_willStay_trials{ne}{n},1),thistrial(sortIdx,:),[minC,maxC])
        colormap(flipud(red2blue))
        title('stay trials')
        hold on
        thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==1);
        hold on;
        plot(thisTimestamp(sortIdx)+3,1:size(out.dff_nrew_willStay_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
        smallcolorbar()
        
        
        
        figure();
        subplot(2,1,1)
        thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==0);
        [~,sortIdx] = sort(thisTimestamp);
        imagesc(xaxis,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),out.dff_nrew_willSwitch_trials{ne}{n}(sortIdx,:),[minC,maxC])
        colormap(flipud(red2blue))
        smallcolorbar()
        hold on;
        plot(thisTimestamp(sortIdx)+3,1:size(out.dff_nrew_willSwitch_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
        title('switch trials')
        subplot(2,1,2)
         thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==1);
        [~,sortIdx] = sort(thisTimestamp);
        imagesc(xaxis,1:size(out.dff_nrew_willStay_trials{ne}{n},1),out.dff_nrew_willStay_trials{ne}{n}(sortIdx,:),[minC,maxC])
        colormap(flipud(red2blue))
        title('stay trials')
        hold on
        thisTimestamp = out.latency{n}(out.reward{n}==0&out.willStay{n}==1);
        hold on;
        plot(thisTimestamp(sortIdx)+3,1:size(out.dff_nrew_willStay_trials{ne}{n},1),'sw','MarkerFaceColor','w','MarkerSize',2)
        smallcolorbar()
        pause()
        close all
        
    end
    
    
  
    
    
end