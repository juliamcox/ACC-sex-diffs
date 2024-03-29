function choiceHist_maleFemale

savehere = whereAreWe('data'); %generate save location
plotParams = load(fullfile(savehere,'plotParams.mat')); 

% load data 
load(fullfile(savehere,'choiceHist_maleFemale.mat'));


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
