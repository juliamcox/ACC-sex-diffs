function [r,p]=opto_control_indVariability(behaviorTable,cohort,ext,zscoreFlag,groupingVar,laserType)

if nargin == 1
   cohort = 'ACC_DMS_nphr';
   ext = 'LAS';
   zscoreFlag = 0; 
   groupingVar = 'none';
   laserType = 1;
end

% Parameters
basefilename    = fullfile(whereAreWe('data'));
% fits: mu (scale), shift, lambda (shape)
plotParams = load(fullfile(basefilename,'plotParams.mat'));


%% Value modulation control v. laser effect
% Load control sessions
aids_f = cat(1,generateAnimalList(sprintf('%s_female',cohort)));
aids_m = cat(1,generateAnimalList(sprintf('%s_male',cohort)));
savehere=fullfile(basefilename,'opto_ctrlSess');
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
load(fullfile(basefilename,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,1,laserType)))


%% Plots 

thisPlot_f_laser = cellfun(@(x) x{2}(1), fits.fits_f);
thisPlot_f  = cellfun(@(x) x{1}(1), fits.fits_f);

thisPlot_m_laser = cellfun(@(x) x{2}(1), fits.fits_m);
thisPlot_m  = cellfun(@(x) x{1}(1), fits.fits_m);

ctrl_f = cellfun(@(x) x{1}(1), fits.fits_f);
ctrl_m = cellfun(@(x) x{1}(1), fits.fits_m);

thisPlot_diff_f = thisPlot_f_laser-thisPlot_f;
thisPlot_diff_m = thisPlot_m_laser - thisPlot_m;

[r.all_individualFit,p.all_individualFit]=corr(cat(2,coef_f(1:numel(thisPlot_diff_f)),coef_m(1:numel(thisPlot_diff_m)))',cat(2,thisPlot_diff_f,thisPlot_diff_m)');
[r.f_individualFit,p.f_individualFit]=corr(coef_f(1:numel(thisPlot_diff_f))',thisPlot_diff_f');
[r.m_individualFit,p.m_individualFit]=corr(coef_m(1:numel(thisPlot_diff_m))',thisPlot_diff_m');


y = cat(1,coef_f',coef_m');
x = cat(1,thisPlot_diff_f',thisPlot_diff_m');
figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
scatterhist(x,y,'Marker','.','MarkerSize',15,...
    'Location','NorthEast','Direction','out','LineWidth',[1,1],'Kernel','on')
g = gca;

figure('Units','inches','Position',[5,5,4.25 2.25]); 
subplot(1,2,1);hold on
boxplot(cat(1,coef_f',coef_m'),'outliersize',8,'symbol','+','notch','off')
ylabel('Value coefficient')
set(gca,'YLim',g.YLim)
subplot(1,2,2)
boxplot(cat(1,thisPlot_diff_f',thisPlot_diff_m'),'outliersize',8,'symbol','+','notch','off')
ylabel('Opto effect')
set(gca,'YLim',g.XLim)
