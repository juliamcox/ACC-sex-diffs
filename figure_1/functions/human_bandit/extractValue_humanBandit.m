function extractValue_humanBandit(fList,cohort,qFileLoc)

% fList: list of all subjects used to fit q-learning model in stan
% qFileLoc: location of stan results

% Load q-model parameters and file list 
load(fullfile(qFileLoc,'sList.mat'));
load(fullfile(qFileLoc,'alphas.mat'));
load(fullfile(qFileLoc,'betas.mat'));
load(fullfile(qFileLoc,'sides.mat'));
load(fullfile(qFileLoc,'stays.mat'));

idx = strfind(qFileLoc, '\');
if isempty(idx)
    idx = strfind(qFileLoc, '/');
end
if isempty(idx)
    idx = 1;
end
savehere = fullfile(whereAreWe('figureCode'),'processed_data','human_bandit',cohort);
if ~isdir(savehere)
    mkdir(savehere)
end
lowCounter = 0;
for ns = 1:numel(fList)
 
    % initialize value structure
    qLearn = struct('QDiff',[],'QRight',[],'QLeft',[],'QTot',[],'probLeft',[],'probRight',[],'choice',[],'reward',[],'QChosenDiff',[],'QChosen',[],'QUnchosen',[]);
    % find subject in sList
    subjIdx   = find(contains(sList,fList{ns}));
    if isempty(subjIdx)
        lowCounter = lowCounter+1;
    else
        % extract parameters for this subject
        thisalpha = alphas(:,subjIdx);
        thisbeta  = betas(:,subjIdx);
        thisside  = sides(:,subjIdx);
        thisstay  = stays(:,subjIdx);
        
        % Load behavior data
        load(fullfile(savehere,fList{ns}));
        numTrials = numel(data.choice);
        % initialize variables for each estimate
        Q_L_est     = zeros(size(alphas,1),numTrials);
        Q_R_est     = zeros(size(alphas,1),numTrials);
        delta       = zeros(size(alphas,1),numTrials); % RPE
        % choice and outcome for this session
        thischoice  = data.choice;
        thisreward  = data.reward;
        
        
        % Estimate Q values
        for nt = 1:numel(thischoice)
            if thischoice(nt) == 0
                delta(:,nt) = thisreward(nt)-Q_R_est(:,nt);
                Q_R_est(:,nt+1) = Q_R_est(:,nt) + thisalpha .* delta(:,nt);
                Q_L_est(:,nt+1) = Q_L_est(:,nt);
            elseif thischoice(nt) == 1
                delta(:,nt) = thisreward(nt)-Q_L_est(:,nt);
                Q_L_est(:,nt+1) = Q_L_est(:,nt) + thisalpha .* delta(:,nt);
                Q_R_est(:,nt+1) = Q_R_est(:,nt);
            elseif isnan(thischoice(nt)) % omitted trial
                Q_L_est(:,nt+1) = Q_L_est(:,nt);
                Q_R_est(:,nt+1) = Q_R_est(:,nt);
            end
        end
        
        % previous choice intercept
        prevChoice           = [0; thischoice(1:end-1) .* 2 - 1];
        intercept_prevChoice = thisstay.*repmat(prevChoice',size(thisstay,1),1);
        % side intercept
        intercept_side       = repmat(thisside,1,numel(thischoice));
        % probability left and right choice
        leftC                = repmat(thisbeta,1,numel(thischoice)).*Q_L_est(:,1:end-1);
        rightC               = repmat(thisbeta,1,numel(thischoice)).*Q_R_est(:,1:end-1);
        prob_left = 1 ./ (1 + exp(rightC - leftC + intercept_prevChoice + intercept_side));
        prob_right = 1 ./ (1 + exp(leftC - rightC + intercept_prevChoice + intercept_side));
        
        qLearn.QDiff            = (mean(Q_R_est(:,1:end-1))-mean(Q_L_est(:,1:end-1)))';
        qLearn.QRight           = mean(Q_R_est(:,1:end-1))';
        qLearn.QLeft            = mean(Q_L_est(:,1:end-1))';
        qLearn.probLeft         = mean(prob_left)';
        qLearn.probRight        = mean(prob_right)';
        QChosen                 = zeros(size(Q_R_est,2),1);
        QChosen(thischoice==1)  = mean(Q_L_est(:,thischoice==1))';
        QChosen(thischoice==0)  = mean(Q_R_est(:,thischoice==0))';
        qLearn.QChosen          = QChosen(1:end-1);
        QUnchosen               = zeros(size(Q_R_est,2),1);
        QUnchosen(thischoice==0)= mean(Q_L_est(:,thischoice==0))';
        QUnchosen(thischoice==1)= mean(Q_R_est(:,thischoice==1))';
        qLearn.QUnchosen        = QUnchosen(1:end-1);
        qLearn.QChosenDiff      = qLearn.QChosen-qLearn.QUnchosen;
        qLearn.QTot             = qLearn.QChosen+qLearn.QUnchosen;
        
        qLearn.choice           = thischoice;
        qLearn.reward           = thisreward;
        
        qLearn.alpha            = thisalpha;
        qLearn.beta            = thisbeta;
        qLearn.side            = thisside;
        qLearn.stay            = thisstay;

        save(fullfile(savehere,sprintf('qLearn_%s',fList{ns})),'qLearn')
    end
end

fprintf('# of subject < 100 trials: %d\n', lowCounter);