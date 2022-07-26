function qValue_estrous

% Generate qLearn structure with value estimates for each session 
% saves qLearn_session_all.mat and qLearn_animal_all.mat in each animal's "Operant" folder

% Data generated by dataForQ_estrous.m
% QDiff is right - left 
% QChosenDiff is chosen - unchosen

%% Parameters
fbasename = fullfile(whereAreWe('figurecode'), 'processed_data', 'dataForQ','estrous_session');

%% load data
load(fullfile(fbasename,'estrous_toStan.mat'))
load(fullfile(fbasename,'estrous_data_bySession.mat'),'flist_all')
stanParams = readtable(fullfile(fbasename,'summary_int_pc.csv'));

aids = cellfun(@(x) x(1).name(2:4),flist_all,'UniformOutput',false);
for na = 1:numel(aids)
 
    fprintf('running %s \n',aids{na});
    filename   = fullfile(whereAreWe('behavior'), aids{na}); % where to save the qlearning results
    probCounter = 0; 
    sessionIDs = find(subj_idx==na);
    thisFlist = {flist_all{na}.name};
   
    qLearn = struct('QDiff',[],'QRight',[],'QLeft',[],'probLeft',[],'probRight',[],'choice',[],'reward',[],'session',[],'trial',[],'QChosenDiff',[],'QChosen',[],'QUnchosen',[], 'estChoice',[]);
    qLearn_mouse = struct('QDiff',[],'QRight',[],'QLeft',[],'probLeft',[],'probRight',[],'choice',[],'reward',[],'session',[],'trial',[],'QChosenDiff',[],'QChosen',[],'QUnchosen',[], 'estChoice', []);
    alpha_mouse = stanParams.mean(strcmp(stanParams.Var1,sprintf('alpha_mice_m[%d]',na)));
    beta_mouse = stanParams.mean(strcmp(stanParams.Var1,sprintf('beta_mice_m[%d]',na)));
    side_mouse = stanParams.mean(strcmp(stanParams.Var1,sprintf('side_mice_m[%d]',na)));
    stay_mouse = stanParams.mean(strcmp(stanParams.Var1,sprintf('stay_mice_m[%d]',na)));
    fList_remove = [];
    for ns = 1:numel(sessionIDs)
        if exist(fullfile(filename,sprintf('%s_%s.mat',thisFlist{ns}(2:end-4),thisFlist{ns}(end-2:end)))) == 0
            fList_remove = cat(1,fList_remove,ns);
        else
        load(fullfile(filename,sprintf('%s_%s.mat',thisFlist{ns}(2:end-4),thisFlist{ns}(end-2:end))));
        
        numTrials = numel(data.choice); 
        %% session level estimates
try        
        alpha(ns) = stanParams.mean(strcmp(stanParams.Var1,sprintf('alphas[%d]',sessionIDs(ns))));
        beta(ns)  = stanParams.mean(strcmp(stanParams.Var1,sprintf('betas[%d]',sessionIDs(ns))));
        stay(ns)  = stanParams.mean(strcmp(stanParams.Var1,sprintf('stays[%d]',sessionIDs(ns))));
        side(ns)  = stanParams.mean(strcmp(stanParams.Var1,sprintf('sides[%d]',sessionIDs(ns))));
catch
    keyboard
end
        if sum(data.choice(data.choice~=-1)' == squeeze(c(na,ns,1:sum(data.choice~=-1)))) ~= sum(data.choice~=-1) || sum(data.reward(data.choice~=-1)' == squeeze((r(na,ns,1:sum(data.choice~=-1))))) ~= sum(data.choice~=-1)
            probCounter = probCounter+1; 
        end
        
        
        
        thisalpha = alpha(ns);
        thisbeta  = beta(ns); 
        thisstay  = stay(ns);
        thisside  = side(ns); 
        % initialize variables for each estimate
        Q_L_est     = zeros(1,numTrials);
        Q_R_est     = zeros(1,numTrials);
        delta       = zeros(1,numTrials); % RPE
        
        % choice and outcome for this session
        thischoice = data.choice';
        thisreward = data.reward';
        
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
            elseif thischoice(nt) == -1
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
        
        estChoice               = arrayfun(@(x) rand(1) < x, mean(prob_right,1));
        
        qLearn.QDiff            = cat(1,qLearn.QDiff, ((Q_R_est(:,1:end-1))-mean(Q_L_est(:,1:end-1)))');
        qLearn.QRight           = cat(1,qLearn.QRight, (Q_R_est(:,1:end-1))');
        qLearn.QLeft            = cat(1,qLearn.QLeft, (Q_L_est(:,1:end-1))');
        qLearn.probLeft         = cat(1,qLearn.probLeft, (prob_left)');
        qLearn.probRight        = cat(1,qLearn.probRight, (prob_right)');
        QChosen                 = zeros(size(Q_R_est,2),1);
        QChosen(thischoice==1)  = (Q_L_est(:,thischoice==1))';
        QChosen(thischoice==0)  = (Q_R_est(:,thischoice==0))'; 
        qLearn.QChosen          = cat(1,qLearn.QChosen,QChosen(1:end-1)); 
        QUnchosen               = zeros(size(Q_R_est,2),1);
        QUnchosen(thischoice==0)= (Q_L_est(:,thischoice==0))';
        QUnchosen(thischoice==1)= (Q_R_est(:,thischoice==1))'; 
        qLearn.QUnchosen        = cat(1,qLearn.QUnchosen,QUnchosen(1:end-1));
        qLearn.choice           = cat(1,qLearn.choice, thischoice);
        qLearn.reward           = cat(1,qLearn.reward,thisreward);
        qLearn.session          = cat(1,qLearn.session, ones(size(thisreward)).*ns);
        qLearn.trial            = cat(1,qLearn.trial, (1:numel(thisreward))');
        qLearn.estChoice        = cat(2,qLearn.estChoice, estChoice); 
        %% mouse level estimates 
%         % initialize variables for each estimate
%         Q_L_est     = zeros(size(alphas,1),numTrials);
%         Q_R_est     = zeros(size(alphas,1),numTrials);
%         delta       = zeros(size(alphas,1),numTrials); % RPE
%        
%         thisalpha = (alphas_mice(:,na));
%         thisbeta  = betas_mice(:,na); 
%         thisstay  = stays_mice(:,na);
%         thisside  = sides_mice(:,na); 
%         
%         % Estimate Q values
%         for nt = 1:numel(thischoice)
%             if thischoice(nt) == 0
%                 delta(:,nt) = thisreward(nt)-Q_R_est(:,nt);
%                 Q_R_est(:,nt+1) = Q_R_est(:,nt) + thisalpha .* delta(:,nt);
%                 Q_L_est(:,nt+1) = Q_L_est(:,nt);
%             elseif thischoice(nt) == 1
%                 delta(:,nt) = thisreward(nt)-Q_L_est(:,nt);
%                 Q_L_est(:,nt+1) = Q_L_est(:,nt) + thisalpha .* delta(:,nt);
%                 Q_R_est(:,nt+1) = Q_R_est(:,nt);
%             elseif thischoice(nt) == -1
%                 Q_L_est(:,nt+1) = Q_L_est(:,nt);
%                 Q_R_est(:,nt+1) = Q_R_est(:,nt);
%             end
%         end
%         
%         % previous choice intercept
%         prevChoice           = [0; thischoice(1:end-1) .* 2 - 1];
%         intercept_prevChoice = thisstay.*repmat(prevChoice',size(thisstay,1),1);
%         % side intercept 
%         intercept_side       = repmat(thisside,1,numel(thischoice));
%         % probability left and right choice
%         leftC                = repmat(thisbeta,1,numel(thischoice)).*Q_L_est(:,1:end-1);
%         rightC               = repmat(thisbeta,1,numel(thischoice)).*Q_R_est(:,1:end-1);
%         prob_left = 1 ./ (1 + exp(rightC - leftC + intercept_prevChoice + intercept_side));
%         prob_right = 1 ./ (1 + exp(leftC - rightC + intercept_prevChoice + intercept_side));
%         estChoice            = arrayfun(@(x) rand(1) < x, mean(prob_right));
% 
%         qLearn_mouse.QDiff            = cat(1,qLearn_mouse.QDiff, (mean(Q_R_est(:,1:end-1))-mean(Q_L_est(:,1:end-1)))');
%         qLearn_mouse.QRight           = cat(1,qLearn_mouse.QRight, mean(Q_R_est(:,1:end-1))');
%         qLearn_mouse.QLeft            = cat(1,qLearn_mouse.QLeft, mean(Q_L_est(:,1:end-1))');
%         qLearn_mouse.probLeft         = cat(1,qLearn_mouse.probLeft, mean(prob_left)');
%         qLearn_mouse.probRight        = cat(1,qLearn_mouse.probRight, mean(prob_right)');
%         QChosen                       = zeros(size(Q_R_est,2),1);
%         QChosen(thischoice==1)        = mean(Q_L_est(:,thischoice==1))';
%         QChosen(thischoice==0)        = mean(Q_R_est(:,thischoice==0))'; 
%         qLearn_mouse.QChosen          = cat(1,qLearn_mouse.QChosen,QChosen(1:end-1)); 
%         QUnchosen                     = zeros(size(Q_R_est,2),1);
%         QUnchosen(thischoice==0)      = mean(Q_L_est(:,thischoice==0))';
%         QUnchosen(thischoice==1)      = mean(Q_R_est(:,thischoice==1))'; 
%         qLearn_mouse.QUnchosen        = cat(1,qLearn_mouse.QUnchosen,QUnchosen(1:end-1));
%         qLearn_mouse.choice           = cat(1,qLearn_mouse.choice, thischoice);
%         qLearn_mouse.reward           = cat(1,qLearn_mouse.reward,thisreward);
%         qLearn_mouse.session          = cat(1,qLearn_mouse.session, ones(size(thisreward)).*ns);
%         qLearn_mouse.trial            = cat(1,qLearn_mouse.trial, (1:numel(thisreward))');
%         qLearn_mouse.estChoice       = cat(2,qLearn_mouse.estChoice, estChoice); 
    end
    end
%     qLearn.beta_mouse = mean(betas_mice(:,na));
%     qLearn.alpha_mouse = mean(alphas_mice(:,na));
%     qLearn.stay_mouse = mean(stays_mice(:,na));
%     qLearn.bias_mouse = mean(sides_mice(:,na));
    qLearn.beta_mouse = beta_mouse;
    qLearn.alpha_mouse = alpha_mouse;
    qLearn.stay_mouse = stay_mouse;
    qLearn.bias_mouse = side_mouse;
    qLearn.alpha_session = alpha;
    qLearn.beta_session = beta;
    qLearn.stay_session = stay;
    qLearn.side_session = side;
    qLearn.QChosenDiff = nan(size(qLearn.QRight));
    qLearn.QChosenDiff(qLearn.choice==0) = qLearn.QRight(qLearn.choice==0)-qLearn.QLeft(qLearn.choice==0);
    qLearn.QChosenDiff(qLearn.choice==1) = qLearn.QLeft(qLearn.choice==1)-qLearn.QRight(qLearn.choice==1);
    qLearn.fList_remove = fList_remove;
    qLearn.fList = flist_all{na};
    
    
    qLearn_mouse.fList = flist_all{na};
    qLearn.probCounter = probCounter;
    
    save(fullfile(filename, 'qLearn_session_summary_2022.mat'),'qLearn');
    qLearn = qLearn_mouse;
  %  save(fullfile(filename, 'qLearn_animal_summary_2022.mat'),'qLearn'); 
    clear beta alpha stay side qLearn qLearn_mouse 
    fprintf('%s prob count: %d \n',aids{na}, probCounter);
end


end

function out = inv_logit(arr)
    %% Elementwise inverse logit (logistic) function
    out = 1 / (1 + exp(-arr));
end

function out = phi_approx(arr)
    %%Elementwise fast approximation of the cumulative unit normal. 
    %%For details, see Bowling et al. (2009). "A logistic approximation 
    %%to the cumulative normal distribution."'''
    out = inv_logit(0.07056 * arr^3 + 1.5976 * arr);
end
