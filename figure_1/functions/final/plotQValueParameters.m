function plotQValueParameters(params_f,params_m)

%%% Figure 1e

% Plot q-value model parameter estimates for males and females (control sessions)

%% Parameters

plotParams    = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat')); % load plot parameters 


%% Plot parameters for males and females


figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
female = cat(1,ones(size(params_f.alpha')),2.*ones(size(params_m.alpha')),3.*ones(size(params_f.alpha')),4.*ones(size(params_m.alpha')),5.*ones(size(params_f.alpha')),6.*ones(size(params_m.alpha')),7.*ones(size(params_f.alpha')),8.*ones(size(params_m.alpha')));
b=boxplot(cat(1,params_f.alpha',params_m.alpha',params_f.beta',params_m.beta',params_f.stay',params_m.stay',params_f.side',params_m.side'),...
    female,'Labels',{'\alpha';'{alpha}';'beta';'beta';'stay';'stay';'side';'side'},'Color',cat(1,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC,plotParams.femaleC,plotParams.maleC),'OutlierSize',8,'Symbol','.','Notch','off');
ax = gca;
ax.TickLabelInterpreter = 'tex';

ax.XTickLabel= {'\alpha';[];'\beta_{value}';[];'\beta_{stay}';[];'\beta_{side}';[]};

box off



%% Make table 
tbl_f = table([1:numel(params_f.alpha)]',repmat({'female'},numel(params_f.alpha),1),'VariableNames',{'Animal';'Sex'});
tbl_f.alpha = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.alpha,params_f.alpha_sem,'UniformOutput',false)';
tbl_f.beta_value = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.beta,params_f.beta_sem,'UniformOutput',false)';
tbl_f.beta_stay = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.stay,params_f.stay_sem,'UniformOutput',false)';
tbl_f.beta_side = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_f.side,params_f.side_sem,'UniformOutput',false)';

tbl_m = table([numel(params_f.alpha)+1:numel(params_m.alpha)+numel(params_f.alpha)]',repmat({'male'},numel(params_m.alpha),1),'VariableNames',{'Animal';'Sex'});
tbl_m.alpha = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.3f')), params_m.alpha,params_m.alpha_sem,'UniformOutput',false)';
tbl_m.beta_value = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.beta,params_m.beta_sem,'UniformOutput',false)';
tbl_m.beta_stay = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.stay,params_m.stay_sem,'UniformOutput',false)';
tbl_m.beta_side = arrayfun(@(x,y) cat(2,num2str(x,'%0.2f'),char(177),num2str(y,'%0.2f')), params_m.side,params_m.side_sem,'UniformOutput',false)';

writetable(tbl_f,fullfile(whereAreWe('figurecode'),'processed_data','stats_tables','qmodel_f.csv'))
writetable(tbl_m,fullfile(whereAreWe('figurecode'),'processed_data','stats_tables','qmodel_m.csv'))

end


