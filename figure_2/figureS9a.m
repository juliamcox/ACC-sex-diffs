function [r,p]=figureS9a(behaviorTable,cohort,ext,zscoreFlag,groupingVar,laserType)

if nargin == 1
   cohort = 'ACC_DMS_nphr';
   ext = 'LAS';
   zscoreFlag = 0; 
   groupingVar = 'none';
   laserType = 1;
end

% Parameters
basefilename    = fullfile(whereAreWe('figurecode'), 'processed_data');
% fits: mu (scale), shift, lambda (shape)
plotParams = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));


%% Value modulation control v. laser effect
% Load control sessions
aids_f = cat(1,generateAnimalList(sprintf('%s_female',cohort)));
aids_m = cat(1,generateAnimalList(sprintf('%s_male',cohort)));
behaviorTable = behaviorTable(boolean(sum(contains(behaviorTable.aID,aids_f)|contains(behaviorTable.aID,aids_m),2)),:);


% effects of value fit separately
for na = 1:numel(aids_m)
    X = behaviorTable(strcmp(behaviorTable.aID, aids_m{na}),:); 
    X = X(X.laserSession==0,:); 
    X.trial = zscore(X.trial); 
    X.Value = (X.qChosenDiff); 
    f = 'trialInit_thresh ~ Value + (1+Value|session)';
    le = fitlme(X,f,'DummyVarCoding','effects');
    coef_m(na) = le.Coefficients.Estimate(contains(le.Coefficients.Name,'Value'));  
end
for na = 1:numel(aids_f)
    X = behaviorTable(strcmp(behaviorTable.aID, aids_f{na}),:); 
    X = X(X.laserSession==0,:); 
    X.trial = zscore(X.trial); 
    X.Value = (X.qChosenDiff); 
    f = 'trialInit_thresh ~ Value + (1+Value|session)';
    le = fitlme(X,f,'DummyVarCoding','effects');
    coef_f(na) = le.Coefficients.Estimate(contains(le.Coefficients.Name,'Value'));  
end


%% Inverse Gaussian 
% Load inverse gaussian fits of opto data
load(fullfile(basefilename,'fig2',sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,1,laserType)))


%% Plots 

thisPlot_f_laser = cellfun(@(x) x{2}(1), fits.fits_f);
thisPlot_f  = cellfun(@(x) x{1}(1), fits.fits_f);

thisPlot_m_laser = cellfun(@(x) x{2}(1), fits.fits_m);
thisPlot_m  = cellfun(@(x) x{1}(1), fits.fits_m);

ctrl_f = cellfun(@(x) x{1}(1), fits.fits_f);
ctrl_m = cellfun(@(x) x{1}(1), fits.fits_m);




figure(); hold on
plot([-4 1],[0 0],'--','Color',[.3 .3 .3],'LineWidth',1.5)
plot([0 0],[3,-4],'--','Color',[.3 .3 .3],'LineWidth',1.5)

thisPlot_diff_f = thisPlot_f_laser-thisPlot_f;
thisPlot_diff_m = thisPlot_m_laser - thisPlot_m;
scatter(thisPlot_diff_f,coef_f(1:numel(thisPlot_diff_f)),30,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none')
scatter(thisPlot_diff_m,coef_m(1:numel(thisPlot_diff_m)),30,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none')
ylabel('Value coefficient')
xlabel('Scale parameter (laser - no laser)')
title('Individual fits')
[r.all_individualFit,p.all_individualFit]=corr(cat(2,coef_f(1:numel(thisPlot_diff_f)),coef_m(1:numel(thisPlot_diff_m)))',cat(2,thisPlot_diff_f,thisPlot_diff_m)');
[r.f_individualFit,p.f_individualFit]=corr(coef_f(1:numel(thisPlot_diff_f))',thisPlot_diff_f');
[r.m_individualFit,p.m_individualFit]=corr(coef_m(1:numel(thisPlot_diff_m))',thisPlot_diff_m');

female = cat(2,zeros(size(thisPlot_diff_f)), ones(size(thisPlot_diff_m)));
y = cat(1,coef_f',coef_m');
x = cat(1,thisPlot_diff_f',thisPlot_diff_m');
figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
scatterhist(x,y,'Group',female,'Color',cat(1,plotParams.femaleC,plotParams.maleC),'Marker','.','MarkerSize',15,...
    'Location','NorthEast','Direction','out','LineWidth',[1,1],'Kernel','on')


