function [fit_m,fit_f] = humanBandit_ageXlatency_scatter(latency_f,latency_m)

plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
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
xlabel('Age (years)')
ylabel('Trial initiation latency (s)'); 

