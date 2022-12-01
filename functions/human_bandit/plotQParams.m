function plotQParams(alpha,alpha_m,beta,beta_m,stay,stay_m,side,side_m)
plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));



figure('Units','inches','Position',[5,5,4.25 2.25]); hold on
female = cat(1,ones(size(alpha')),2.*ones(size(alpha_m')));
pOrig=subplot(1,4,1);
b=boxplot(cat(1,alpha',alpha_m'),female,'Color',cat(1,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
ax = gca;
ax.Position(3:4) = [0.1502 0.7134];
%ax.TickLabelInterpreter = 'tex';
xlabel({'\alpha';[]});
ax.XTick = [];
box off
set(gca,'YLim',[0 1])

p=subplot(1,4,2); p.Position(3) = pOrig.Position(3);
b=boxplot(cat(1,beta',beta_m'),female,'Color',cat(1,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
ax = gca;
ax.Position(3:4) = [0.1502 0.7134];
xlabel({'\beta_{value}';[]});
ax.YLim(1) = 0;
ax.XTick = [];
box off
ax.Position(3:4) = [0.1502 0.7134];

p=subplot(1,4,3);
p.Position(3) = pOrig.Position(3);
b=boxplot(cat(1,stay',stay_m'),female,'Color',cat(1,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
ax = gca;
ax.Position(3:4) = [0.1502 0.7134];
xlabel({'\beta_{stay}';[]});
ax.XTick = [];box off


p=subplot(1,4,4); p.Position(3) = pOrig.Position(3);
b=boxplot(cat(1,side',side_m'),female,'Color',cat(1,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
ax = gca;
ax.Position(3:4) = [0.1502 0.7134];
ax.XTick = [];
xlabel({'\beta_{side}';[]});
box off
