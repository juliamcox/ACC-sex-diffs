function [fit_m,fit_f] = humanBandit_ageXlatency_scatter(latency_f,latency_m)

plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
figure(), scatter(latency_f.age,cellfun(@(x) nanmean(x),latency_f.trialStart_trials),20,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.50)
hold on
scatter(latency_m.age,cellfun(@(x) nanmean(x),latency_m.trialStart_trials),20,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.5)

f = 'latency ~ age';
T = table(latency_f.age','VariableNames',{'age'});
T.latency = cellfun(@(x) nanmean(x),latency_f.trialStart_trials)';
fit_f = fitlm(T,f);

f = 'latency ~ age';
T = table(latency_m.age','VariableNames',{'age'});
T.latency = cellfun(@(x) nanmean(x),latency_m.trialStart_trials)';
fit_m = fitlm(T,f);

intercept = fit_f.Coefficients.Estimate(1);
age = fit_f.Coefficients.Estimate(2);
plot([19:71], arrayfun(@(x) intercept + age*x, [19:71]), 'Color',plotParams.femaleC,'LineWidth',1.5)
intercept = fit_m.Coefficients.Estimate(1);
age = fit_m.Coefficients.Estimate(2);
plot([19:71], arrayfun(@(x) intercept + age*x, [19:71]), 'Color',plotParams.maleC,'LineWidth',1.5)


x = cat(1,latency_m.age',latency_f.age');
y = cat(1,cellfun(@(x) nanmean(x),latency_m.trialStart_trials)',cellfun(@(x) nanmean(x),latency_f.trialStart_trials)');
female = cat(1,zeros(size(latency_m.age')),ones(size(latency_f.age')));
figure();
scatterhist(x,y,'Group',female,'Color',cat(1,plotParams.maleC,plotParams.femaleC),'Marker','.','MarkerSize',15,...
'Location','NorthEast','Direction','out','LineWidth',[1,1],'Kernel','on')