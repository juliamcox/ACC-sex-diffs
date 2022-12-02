function lme = pairedEPSP(data) 


plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
%% Plot 
tbl_m = data;



pairs = unique(tbl_m.pair);
for np = 1:numel(pairs)
   d1_m(np) = tbl_m.current(tbl_m.pair==pairs(np)&tbl_m.msn_type==1); 
   d2_m(np) = tbl_m.current(tbl_m.pair==pairs(np)&tbl_m.msn_type==2);
end


figure()

subplot(1,2,2)
plot([d1_m;d2_m],'Color',plotParams.maleC,'LineWidth',1)
box off
set(gca,'XLim', [.5 2.5], 'YLim', [0 25], 'XTick',[1 2],'XTickLabel',{'D1R';'D2R'})
ylabel('EPSP amplitude (mV)')
%% LME
f = 'current ~ msn_type + (1|pair:mouse) + (1|mouse)';

data.msn_type = categorical(data.msn_type);
data.pair = categorical(data.pair);
data.mouse = categorical(data.mouse); 

lme = fitlme(data,f,'DummyVarCoding','effect');