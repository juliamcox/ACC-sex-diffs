Sfunction figure1F_plotExampleSession(aids)

%%% Generate example behavior session with estimated Q-values and choice for figure 1F 

% Required functions: qValue_session.m (by session fit) or qLearning_transformation
    % (by animal fit); TABGUI for generation of behavior files 

% aids: list of animal IDs

qFile = 'qLearn_session_all_2022.mat'; % which value file 
fbasename = fullfile(whereAreWe('behavior')); 

for na = 1:numel(aids)
    close all
    % Load value estimates
    load(fullfile(fbasename,aids{na},qFile));
        for nf = 1:numel(qLearn.fList)
            if contains(qLearn.fList(nf).name, 'TAB') % only plot examples from control sessions
                %%% Load behavior file
                idx = strfind(qLearn.fList(nf).name,'.');
                thisFilename = sprintf('%s_%s.mat',qLearn.fList(nf).name(2:idx-1),qLearn.fList(nf).name(idx+1:end));
                load(fullfile(fbasename,aids{na},thisFilename));
                %%% Choose 100 trials to plot 
                trialIdx                        = randi(numel(data.choice)-150,1);
                trialIdx                        = trialIdx:trialIdx+149;
                %%% Extract variables for plotting
                % estimated relative value (right - left) 
                qDiff                                = nan(size(data.choice));
                qDiff(data.choice~=-1)               = qLearn.QRight(qLearn.session==nf) - qLearn.QLeft(qLearn.session==nf);
                qDiff                                = qDiff(trialIdx); 
                % predicted choice
                predChoice                           = nan(size(data.choice));
                probLeft                             = nan(size(data.choice));
                probLeft(data.choice~=-1)            = qLearn.probLeft(qLearn.session==nf);
                probRight                            = nan(size(data.choice)); 
                probRight(data.choice~=-1)           = qLearn.probRight(qLearn.session==nf);
                predChoice(rand(1,numel(predChoice)) < probRight) = 0;
                predChoice(isnan(predChoice))        = 1;
                predChoice(probLeft<probRight)       = 0;
                predChoice                           = predChoice(trialIdx); % switch right choice to -1 
                predChoice(predChoice==0)            = -1;
                % mouse's choice
                choice                               = data.choice(trialIdx);
                choice(choice==-1)                   = NaN;
                choice                               = choice.*2 - 1; % switch right choice to -1 
                % reward
                reward                               = nan(size(data.reward));
                reward(data.choice==0&data.reward==1)= -1.1; % right reward -1 
                reward(data.choice==1&data.reward==1)= 1.1; 
                reward                               = reward(trialIdx); 
                % high probability block
                blockID                              = data.blocks.blockIDAll(trialIdx);
                rightBlock                           = nan(size(blockID));
                leftBlock                            = nan(size(blockID));
                rightBlock(blockID==1)               = -1.75;
                leftBlock(blockID==2)                = 1.75;
                probRight                            = probRight(trialIdx); 
                
                figure(); hold on
                xaxis = 1:numel(trialIdx);
                plot(xaxis, reward, '.','MarkerSize',8)
                plot(xaxis, choice.*1.25, 's','MarkerSize',3, 'MarkerFaceColor','k','MarkerEdgeColor','none')
                %plot(xaxis, predChoice.*1.5, 's','MarkerSize',3, 'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor','none')
                plot(xaxis, probRight,'LineWidth',2)
                plot(xaxis, qDiff, 'Color', [255 54 7]./255)
                plot(xaxis, rightBlock,'LineWidth',2)
                plot(xaxis, leftBlock,'LineWidth',2)
                plot(1:50, ones(1,50).*-2.1, 'k', 'LineWidth',1) %trial scale 
                set(gca,'YLim', [-2 2])
                title(sprintf('%s session %s Trials %s:%s',aids{na}, num2str(nf), num2str(trialIdx(1)),num2str(trialIdx(end))))
                if sum(isnan(choice)) ~= 0
                    fprintf('there are omitted trials in this session \n')
                end
                pause()
                clear trialIdx qDiff predChoice choice reward blockID
            end
        end
  
    
end

    