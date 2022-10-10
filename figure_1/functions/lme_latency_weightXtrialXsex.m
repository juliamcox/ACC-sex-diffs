function stats = lme_latency_weightXtrialXsex(aids_f,aids_m,behaviorTable,valType)

for nv = 1:numel(valType)
%% Mixed-effects model with weight, value, trial number and sex 
% select animals
idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
X = behaviorTable(cell2mat(idx),:);
% select non-laser sessions
X = X(X.laserSession==0,:);
X.trial = nanzscore(X.trial); 
X.female = categorical(X.female);
eval(sprintf('X.%s_quant = categorical(X.%s_quant);',valType{nv},valType{nv}));

f = sprintf('trialInit_thresh ~ %s_quant*female*weight_zscore + %s_quant*female*trial + (1+%s_quant*trial + %s_quant*weight_zscore|aID)+(1+%s_quant*trial|aID:session)',valType{nv},valType{nv},valType{nv},valType{nv},valType{nv});

mdl = fitlme(X,f,'DummyVarCoding','effects');
stats.mdl = mdl;
stats.mdl_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite')); 
stats.mdl_coeff = dataset2cell(mdl.Coefficients);
weights = min(X.weight_zscore):.1:max(X.weight_zscore);
trial = unique(X.trial); 
% posthoc contrasts (male v. female for each value bin); mid session, across weights
for nb = 1:max(behaviorTable.qChosenDiff_quant)-1
   h0 = zeros(1,size(stats.mdl_coeff,1)-1);
   thisIdx = find(contains(stats.mdl_coeff(:,1), 'female')&contains(stats.mdl_coeff(:,1),num2str(nb)));
   h0(thisIdx-1) = 2; 
   h0(3) = 2; 
   thisIdx = find(contains(stats.mdl_coeff(:,1),'trial'));
   h0(thisIdx-1) =  h0(thisIdx-1).*nanmean(X.trial); 
   h = [];
   for nw = 1:numel(weights)
       thisIdx = find(contains(stats.mdl_coeff(:,1),'weight'));
       temp = h0;
       temp(thisIdx-1) = temp(thisIdx-1).*weights(nw);
       h = cat(1, h, temp);
   end  
   [stats.mdl_coefTest_p(nb),stats.mdl_coefTest_f(nb),stats.mdl_coefTest_df1(nb),stats.mdl_coefTest_df2(nb)] = coefTest(stats.mdl,h);
end
h0 = zeros(1,size(stats.mdl_coeff,1)-1);
thisIdx = find(contains(stats.mdl_coeff(:,1), 'female')&contains(stats.mdl_coeff(:,1),'qChosenDiff_quant'))-1;
h0(thisIdx) = -2;
h0(3) = 2;
thisIdx = find(contains(stats.mdl_coeff(:,1),'trial'));
h0(thisIdx-1) =  h0(thisIdx-1).*nanmean(X.trial);
for nw = 1:numel(weights)
    thisIdx = find(contains(stats.mdl_coeff(:,1),'weight'));
    temp = h0;
    temp(thisIdx-1) = temp(thisIdx-1).*weights(nw);
    h = cat(1, h, temp);
end
[stats.mdl_coefTest_p(4),stats.mdl_coefTest_f(4),stats.mdl_coefTest_df1(4),stats.mdl_coefTest_df2(4)] = coefTest(stats.mdl,h);

end

%% Mixed effects regression for outcome 

% % select animals
% idx = cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_f,'UniformOutput',false);
% idx = cat(1,idx,cellfun(@(x) find(strcmp(behaviorTable.aID,x)),aids_m,'UniformOutput',false));
% X = behaviorTable(cell2mat(idx),:);
% % select non-laser sessions
% X = X(X.laserSession==0,:);
% X.trial = nanzscore(X.trial); 
% X.female = categorical(X.female);
% X.previousReward = categorical(X.previousReward); 
% X.aID = categorical(X.aID);
% X.session = categorical(X.session);
% 
% 
% f = 'trialInit ~ previousReward*female  + (1+ previousReward|aID)';
% 
% mdl = fitlme(X,f,'DummyVarCoding','effects');
% stats.mdlOut = mdl;
% stats.mdlOut_anova = dataset2cell(anova(mdl,'DFMethod','Satterthwaite')); 
% stats.mdlOut_coeff = dataset2cell(mdl.Coefficients);

% h_m_rew = [1 -1 1 -1];
% h_m_nrew = [1 1 1 1];
% h_f_rew = [1 -1 -1 1];
% h_f_nrew = [1 1 -1 -1];
% h = h_m_rew-h_f_rew;
% [stats.contrasts.rew_p_p,stats.contrasts.rew_p_f,stats.contrasts.rew_p_df1,stats.contrasts.rew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_m_nrew-h_f_nrew;
% [stats.contrasts.nrew_p_p,stats.contrasts.nrew_p_f,stats.contrasts.nrew_p_df1,stats.contrasts.nrew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_f_rew-h_f_nrew;
% [stats.contrasts.f_p_p,stats.contrasts.f_p_f,stats.contrasts.f_p_df1,stats.contrasts.f_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_m_rew-h_m_nrew;
% [stats.contrasts.m_p_p,stats.contrasts.m_p_f,stats.contrasts.m_p_df1,stats.contrasts.m_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');


% h_m_rew = [1 -1 0 1 0 0 -1 0 0 0 0 0];
% h_m_nrew = [1 1 0 1 0 0 1 0 0 0 0 0];
% h_f_rew = [1 -1 0 -1 0 0 1 0 0 0 0 0];
% h_f_nrew = [1 1 0 -1 0 0 -1 0 0 0 0 0];
% h = h_m_rew-h_f_rew;
% [stats.contrasts.rew_p_p,stats.contrasts.rew_p_f,stats.contrasts.rew_p_df1,stats.contrasts.rew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_m_nrew-h_f_nrew;
% [stats.contrasts.nrew_p_p,stats.contrasts.nrew_p_f,stats.contrasts.nrew_p_df1,stats.contrasts.nrew_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_f_rew-h_f_nrew;
% [stats.contrasts.f_p_p,stats.contrasts.f_p_f,stats.contrasts.f_p_df1,stats.contrasts.f_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% h = h_m_rew-h_m_nrew;
% [stats.contrasts.m_p_p,stats.contrasts.m_p_f,stats.contrasts.m_p_df1,stats.contrasts.m_p_df2]   = coefTest(mdl,h,0,'DFMethod','Satterthwaite');
% 
