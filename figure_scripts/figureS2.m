clear 
close all

%% Parameters
cohort          = 'ACC_DMS_nphr';
ext             = 'LAS'; 
zscoreFlag      = 0; % zscore latencies?
basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data');
groupingVar     = 'none';
ver             = 1; %parameterization 
%% Supp Figure 2
laserType       = 2;

load(fullfile(basefilename,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
stats_invGauss_nosePoke=plotInverseGaussian_parameters(fits,{'f';'m'});
plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});


laserType       = 3;
load(fullfile(basefilename,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)))
stats_invGauss_ITI=plotInverseGaussian_parameters(fits,{'f';'m'});
plotInverseGaussian_parameters(fits,{'f_yfp';'m_yfp'});
