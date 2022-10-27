function outcomeScatter_maleFemale(recs_f,recs_m,frameRate,rawFlag,zscoreFlag,epoch)
if nargin < 8
    expt = 'ACC_DMS_imaging';
end


fbasename = fullfile(whereAreWe('imaging'),expt);


savehere = fullfile(whereAreWe('figurecode'),'processed_data','fig4'); %generate save location
if ~isfolder(savehere)
    mkdir(savehere)
end

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); 

try load(fullfile(savehere,['outcomeHist_maleFemale_' epoch '.mat']));
catch
    % Find p-vals and average outcome activity
    switch epoch
        case 'early'
            %early ITI (histLength)
            histLength = [0 3];
            [rewHists_m,nrewHists_m,rew_avg_m,nrew_avg_m,pvals_all_m] = outcomeActivity(recs_m,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
            [rewHists_f,nrewHists_f,rew_avg_f,nrew_avg_f,pvals_all_f] = outcomeActivity(recs_f,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
        case 'full' %fullITI
            histLength = [0];
            [rewHists_m,nrewHists_m,rew_avg_m,nrew_avg_m,pvals_all_m] = outcomeActivity(recs_m,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
            [rewHists_f,nrewHists_f,rew_avg_f,nrew_avg_f,pvals_all_f] = outcomeActivity(recs_f,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
        case 'late'
            histLength = [1 0];
            [rewHists_m,nrewHists_m,rew_avg_m,nrew_avg_m,pvals_all_m] = outcomeActivity(recs_m,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
            [rewHists_f,nrewHists_f,rew_avg_f,nrew_avg_f,pvals_all_f] = outcomeActivity(recs_f,frameRate,rawFlag,zscoreFlag,histLength,fbasename);
    end
end

save(fullfile(savehere,['outcomeHist_maleFemale_' epoch '.mat']));



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

minD = -1;
maxD = 1; 

hold on
[a,b] = histcounts(diff_avg_male,linspace(minD,maxD,30),'Normalization','probability');
[a2,b2] = histcounts(diff_avg_female,linspace(minD,maxD,30),'Normalization','probability');

plot(b(1:end-1),a,'Color',maleC)
plot(b(1:end-1),a2,'Color',femaleC)
plot(median(diff_avg_male),max(cat(2,a,a2))+.02,'v','Color',maleC,'MarkerSize',6)
plot(median(diff_avg_female),max(cat(2,a,a2))+.02,'v','Color',femaleC,'MarkerSize',6)
[p_all,~,stats] = ranksum(diff_avg_male,diff_avg_female);

if p_all<.01
    title('All active neurons *')
else
    title('All active neurons')
end

xlabel('(Reward - No reward activity)/(Reward + No reward activity)')
ylabel('Proportion of neurons')


end

