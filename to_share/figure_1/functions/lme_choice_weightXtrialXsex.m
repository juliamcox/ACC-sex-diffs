function stats = lme_choice_weightXtrialXsex(aids_f,aids_m,behaviorTable,valType)
behaviorTable = behaviorTable(behaviorTable.laserSession == 0,:);
for nv = 1:numel(valType)
%% Mixed-effects model with weight, value, trial number and sex 
% select animals
idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
X = behaviorTable(cell2mat(idx),:);
% select non-laser sessions
X = X(X.laserSession==0,:);
X = X(~isnan(X.trialInit_thresh),:);
X.trial = nanzscore(X.trial); 
X.female = categorical(X.female);
eval(sprintf('X.%s_quant_choice = categorical(X.%s_quant_choice);',valType{nv},valType{nv}));

f = sprintf('choice ~ %s_quant_choice*female*trial + (1+%s_quant_choice*trial|aID)',valType{nv},valType{nv});

mdl = fitglme(X,f,'DummyVarCoding','effects','Distribution','binomial');
eval(sprintf('stats.mdl_%s_quant = mdl;', valType{nv})); 
eval(sprintf('stats.mdl_anova_%s_quant = dataset2cell(anova(mdl));',valType{nv})); 
eval(sprintf('stats.mdl_coeff_%s_quant = dataset2cell(mdl.Coefficients);',valType{nv})); 


end

%% Mixed effects regression for outcome 

% select animals
idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
X = behaviorTable(cell2mat(idx),:);
% select non-laser sessions
X = X(X.laserSession==0,:);
X = X(~isnan(X.trialInit_thresh),:);
X.trial = nanzscore(X.trial); 
X.female = categorical(X.female);
X.previousReward = categorical(X.previousReward); 

f = 'stay ~ previousReward*female*trial + (1+previousReward*trial|aID)';

mdl = fitglme(X,f,'DummyVarCoding','effects','Distribution','Binomial');
stats.mdlOut = mdl;
stats.mdlOut_anova = dataset2cell(anova(mdl)); 
stats.mdlOut_coeff = dataset2cell(mdl.Coefficients);
