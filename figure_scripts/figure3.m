clear all
close all

%% Figure 3d

load(fullfile(whereAreWe('data'),'chr2_current_sex_comparison_matrix.mat'))

lme = chr2_paired_current_lme(chr2currentsexcomparison);

anovaIn = dataset2cell(anova(lme,'DFMethod','Satterthwaite')); 
coeffIn = dataset2cell(lme.Coefficients);
tbl = statsTable(anovaIn, coeffIn, 1, 0,fullfile(whereAreWe('data'),'stats_tables','ST10.csv'));

