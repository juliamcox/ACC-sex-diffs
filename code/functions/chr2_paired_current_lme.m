function [lme] = chr2_paired_current_lme(chr2currentsexcomparison)
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));

tbl = chr2currentsexcomparison;

%Make categorical variables
tbl.cell = categorical(tbl.cell);
tbl.mouse = categorical(tbl.mouse);
tbl.pair = categorical(tbl.pair);
tbl.sex = categorical(tbl.sex);
tbl.msn_type = categorical(tbl.msn_type);


f_ids = unique(tbl.mouse(tbl.sex == categorical(2)));
m_ids = unique(tbl.mouse(tbl.sex == categorical(1)));


%% Perform linear mixed effects model analysis

lme = fitlme(tbl,'current ~ sex * msn_type + (1|mouse) + (1|pair:mouse)','DummyVarCoding','effects');

%% Plot 
tbl = chr2currentsexcomparison;

tbl_f = tbl(tbl.sex==2,:);
pairs = unique(tbl_f.pair);
for np = 1:numel(pairs)
   d1_f(np) = tbl_f.current(tbl_f.pair==pairs(np)&tbl_f.msn_type==1); 
   d2_f(np) = tbl_f.current(tbl_f.pair==pairs(np)&tbl_f.msn_type==2);
end

tbl_m = tbl(tbl.sex==1,:);
pairs = unique(tbl_m.pair);
for np = 1:numel(pairs)
   d1_m(np) = tbl_m.current(tbl_m.pair==pairs(np)&tbl_m.msn_type==1); 
   d2_m(np) = tbl_m.current(tbl_m.pair==pairs(np)&tbl_m.msn_type==2);
end


figure()
% plot females
subplot(1,2,1)
plot([-d1_f;-d2_f],'Color',plotParams.femaleC,'LineWidth',1)
box off
set(gca,'XLim', [.5 2.5], 'YLim', [0 1000], 'XTick',[1 2],'XTickLabel',{'D1R';'D2R'})
ylabel('EPSC amplitude (pA)')
subplot(1,2,2)
plot([-d1_m;-d2_m],'Color',plotParams.maleC,'LineWidth',1)
box off
set(gca,'XLim', [.5 2.5], 'YLim', [0 1000], 'XTick',[1 2],'XTickLabel',{'D1R';'D2R'})
ylabel('EPSC amplitude (pA)')

