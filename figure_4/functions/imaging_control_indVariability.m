function [r,p]=imaging_control_indVariability(perfThresh,qFile,sessionLength,binNum,binNum_choice,intervalThresh,cohort_opto,recs_f,recs_m)


idx = strfind(recs_f,'\');
if isempty(idx{1})
    idx = strfind(recs_f,'/');
end
aids_f = cellfun(@(x,y) x(1:y-1), recs_f,idx,'UniformOutput',false);
idx = strfind(recs_m,'\');
if isempty(idx{1})
    idx = strfind(recs_m,'/');
end
aids_m = cellfun(@(x,y) x(1:y-1), recs_m,idx,'UniformOutput',false);

% Parameters
basefilename    = fullfile(whereAreWe('figureCode'), 'processed_data');
plotParams      = load(fullfile(whereAreWe('figurecode'),'general_code','plotParams.mat'));
frameRate = 20;
zscoreFlag = 1; 
rawFlag = 1; 
%% Value modulation control v. laser effect
% Load control sessions
savehere=fullfile(basefilename,'imaging_ctrlSess');


load(fullfile(basefilename,sprintf('allBehaviorData_perfThresh%s_bin%s_binChoice%s_cohortOpto_%s_%s_intThresh%s_%s',num2str(perfThresh),num2str(binNum),num2str(binNum_choice),cohort_opto, sessionLength,num2str(intervalThresh),qFile)),'behaviorTable');


% effects of value fit separately
for na = 1:numel(aids_m)
    X = behaviorTable(strcmp(behaviorTable.aID,aids_m{na}),:);
    f = 'trialInit ~ qChosenDiff + (1+qChosenDiff|session)';
    le = fitlme(X,f,'DummyVarCoding','effects');
    coef_m(na) = le.Coefficients.Estimate(contains(le.Coefficients.Name,'qChosenDiff')); 
    [randEff,randEffNames] = randomEffects(le);
    idx = strfind(recs_m{na},'_');
    recDate = recs_m{na}(idx(1)+1:idx(2)-1);
    recDate = datenum(recDate,'yyyymmdd');
    sessionDates = unique(X.sessionDate);
    imagingIdx = find(sessionDates==recDate);
    coef_m_rand(na) = coef_m(na)+randEff(strcmp(randEffNames.Name,'qChosenDiff')&strcmp(randEffNames.Level,num2str(imagingIdx)));
end
for na = 1:numel(aids_f)
    X = behaviorTable(strcmp(behaviorTable.aID,aids_f{na}),:);
    f = 'trialInit ~ qChosenDiff + (1+qChosenDiff|session)';
    le = fitlme(X,f,'DummyVarCoding','effects');
    coef_f(na) = le.Coefficients.Estimate(contains(le.Coefficients.Name,'qChosenDiff')); 
    [randEff,randEffNames] = randomEffects(le);
    idx = strfind(recs_f{na},'_');
    recDate = recs_f{na}(idx(1)+1:idx(2)-1);
    recDate = datenum(recDate,'yyyymmdd');
    sessionDates = unique(X.sessionDate);
    imagingIdx = find(sessionDates==recDate);
    coef_f_rand(na) = coef_f(na)+randEff(strcmp(randEffNames.Name,'qChosenDiff')&strcmp(randEffNames.Level,num2str(imagingIdx)));
end


%% Load imaging data and find percent outcome, reward and no reward by animal 

ver = 'outcome';
whichEvent = 6; %no reward 

a = .01; %significace level

fname =  sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s',frameRate,rawFlag,zscoreFlag,ver);
fbasename = fullfile(whereAreWe('imaging'));
fbasename_bs = fullfile(whereAreWe('figurecode'),'general_code', 'basis_sets');

% load basis set and events
[cons ~, ~, ~, ~, ~, bsIDs,~,~] = getEvents(ver,frameRate);


load(fullfile(fbasename_bs, ['bs_' bsIDs '.mat']))
bs = (full(eval(['bs_' bsIDs])));

b_bs_all = cell(1,numel(cons));
a = 0.01; 
for nr = 1:numel(recs_f)
    % load coefficients 
    load(fullfile(fbasename,recs_f{nr},fname),'pvals','b','con_iden');
    % load index of active neurons 
    load(fullfile(fbasename,recs_f{nr},'activeNeurons.mat'))
    b = cell2mat(b);
    b_female = b(2:end,activeIdx);
    pvals = pvals(activeIdx);
    pvals = cell2mat(pvals');
    
    for nn = 1:size(b_female,2)
        for ne = whichEvent
            thisWeights = b_female(con_iden==ne,nn);
            if iscell(bs)
                tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*bs{ne},2)';
            else
                tempWeights = sum(repmat(thisWeights',size(bs,1),1).*bs,2)';
            end
            b_bs_all{ne} = tempWeights;
            clear tempWeights
        end
        aoc_f(nn,ne) = trapz(b_bs_all{ne}');
    end
    propPos(nr,ne) = mean(aoc_f(:,ne)<0&pvals(:,ne)<a);
    propNeg(nr,ne) = mean(aoc_f(:,ne)>0&pvals(:,ne)<a);
    propOut(nr,ne) = mean(pvals(:,ne)<a); 
    
    numPos(nr,ne) = sum(aoc_f(:,ne)<0&pvals(:,ne)<a);
    numNeg(nr,ne) = sum(aoc_f(:,ne)>0&pvals(:,ne)<a);
    numOut(nr,ne) = sum(pvals(:,ne)<a); 
    tot(nr,ne)  = numel(pvals(:,ne));
    aoc_mean(nr,ne) = mean(aoc_f(pvals(:,ne)<a,ne));
    clear aoc_f
end

for nr = 1:numel(recs_m)
    % load coefficients 
    load(fullfile(fbasename,recs_m{nr},fname),'pvals','b','con_iden');
    % load index of active neurons 
    load(fullfile(fbasename,recs_m{nr},'activeNeurons.mat'))
    b = cell2mat(b);
    b_male = b(2:end,activeIdx);
    pvals = pvals(activeIdx);
    pvals = cell2mat(pvals');
    
    for nn = 1:size(b_male,2)
        for ne = whichEvent
            thisWeights = b_male(con_iden==ne,nn);
            if iscell(bs)
                tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*bs{ne},2)';
            else
                tempWeights = sum(repmat(thisWeights',size(bs,1),1).*bs,2)';
            end
            b_bs_all{ne} = tempWeights;
            clear tempWeights
        end
        aoc_m(nn,ne) = trapz(b_bs_all{ne}');
    end
    propPos_m(nr,ne) = mean(aoc_m(:,ne)<0&pvals(:,ne)<a);
    propNeg_m(nr,ne) = mean(aoc_m(:,ne)>0&pvals(:,ne)<a);
    propOut_m(nr,ne) = mean(pvals(:,ne)<a); 
    
    numPos_m(nr,ne) = sum(aoc_m(:,ne)<0&pvals(:,ne)<a);
    numNeg_m(nr,ne) = sum(aoc_m(:,ne)>0&pvals(:,ne)<a);
    numOut_m(nr,ne) = sum(pvals(:,ne)<a); 
    tot_m(nr,ne)  = numel(pvals(:,ne));
    
    aoc_mean_m(nr,ne) = mean(aoc_m(pvals(:,ne)<a,ne));
    clear aoc_m
end


%% Plots 

y = cat(2,coef_f,coef_m);
x = cat(1,(numNeg(:,whichEvent)-numPos(:,whichEvent))./tot(:,whichEvent),(numNeg_m(:,whichEvent)-numPos_m(:,whichEvent))./tot_m(:,whichEvent));

figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
scatterhist(x,y,...
    'Location','NorthEast','Direction','out','LineWidth',[1,1],'Kernel','on');
g = gca;
ylabel('Value modulation')
xlabel('%no rew - % rew')

figure('Units','inches','Position',[5,5,2.25 2.25]); hold on
subplot(1,2,1);
boxplot(y);
ylabel('value modulation')
set(gca,'YLim', g.YLim)
subplot(1,2,2);
boxplot(x);
ylabel('%no rew - % rew')
set(gca,'YLim', g.XLim)

[r.propDiff,p.propDiff] = corr(x,y');


x(y==max(abs(y))) = [];
y(y==max(abs(y))) = [];
[r.diffProp_noOutlier,p.diffProp_noOutlier] = corr(x,y');


