function fit = bilinearEncodingModel_fit(fit,y)

maxIter   = 10;
tolFun    = 1e-3; 



%% Generate predictors 

for ne = 1:numel(fit.eventNames)
    thisEvent = fit.event_times{ne};
    thisEvent(isnan(thisEvent)) = [];
    kernelPreds{ne} = zeros(size(y));
    kernelPreds{ne}(thisEvent-fit.eventDur(ne,1)*fit.frameRate) = 1;
    for ng = 1:numel(fit.gainNames{ne})+1
        thisEvent = fit.event_times{ne};
        thisGain = fit.gain_preds{ne,ng};
        thisGain(isnan(thisEvent)) = [];
        thisEvent(isnan(thisEvent)) = [];
        gainPreds{ne,ng} = zeros(size(y));
        gainPreds{ne,ng}(thisEvent-fit.eventDur(ne,1)*fit.frameRate) = thisGain;
    end
end

%% Fit model with all data

ni = 1;
fchange = inf;

while (ni <= maxIter)
    %%% step one - fit kernel profiles %%%
    X = [];
    for ne = 1:numel(fit.eventNames)
        clear coeff
        % if not the first iteration, find the trial coefficients * trial predictors for each kernel
        if ni > 1
            for ng = 1:numel(fit.gainNames{ne})+1
                coeff(:,ng) = model_parameters.trial_coefficients{ne}(ng).*fit.gain_preds{ne,ng}(~isnan(fit.event_times{ne}));
            end
            coeff = sum(coeff,2);          
        else
            coeff = 1;
        end
        thisPred = kernelPreds{ne};
        thisPred(thisPred==1) = coeff;
        for nb = 1:fit.nbasis(ne)
            temp = conv(thisPred,full(fit.bs{ne}(:,nb)));
            X = cat(2,X,temp(1:numel(y)));
        end
        
    end
    
    thisY = y(sum(X,2)~=0);
    thisX = X(sum(X,2)~=0,:);
    mdl1 = fitlm(thisX,thisY,'linear','intercept',false); 

    idx = 1;
    for ne = 1:numel(fit.eventNames)
        model_parameters.kernel_profiles{ne} = mdl1.Coefficients.Estimate(idx:idx+fit.nbasis(ne)-1);
        idx = idx+fit.nbasis(ne);
        model_parameters.kernel_profiles{ne} = model_parameters.kernel_profiles{ne}./norm(fit.bs{ne}*model_parameters.kernel_profiles{ne},2); % normalize w/ euclidean norm of the kernel
    end
      
     %%% step two - fit trial-by-trial coefficients %%%
     
     X2 = [];
     for ne = 1:numel(fit.eventNames)
         for ng = 1:numel(fit.gainNames{ne})+1
            temp = conv(gainPreds{ne,ng}, fit.bs{ne}*model_parameters.kernel_profiles{ne});
            X2 = cat(2,X2,temp(1:numel(y))); 
         end
     end
     
     thisY = y(sum(X2,2)~=0);
     thisX = X2(sum(X2,2)~=0,:);
     
     
     mdl2 = fitlm(thisX,thisY,'linear','intercept',false); 
     idx = 1; 
     for ne = 1:numel(fit.eventNames)
         model_parameters.trial_coefficients{ne} = mdl2.Coefficients.Estimate(idx:idx+numel(fit.gainNames{ne})+1-1);
         idx = idx+numel(fit.gainNames{ne})+1;
     end
     
     %%% step three - cross-validate %%%
     X = [];
     
     for ne = 1:numel(fit.eventNames)
         clear coeff
         for ng = 1:numel(fit.gainNames{ne})+1
             coeff(:,ng) = model_parameters.trial_coefficients{ne}(ng).*fit.gain_preds{ne,ng}(~isnan(fit.event_times{ne}));
         end
         coeff = sum(coeff,2);   
         
         thisPred = kernelPreds{ne};
         thisPred(thisPred==1) = coeff; 
         
          for nb = 1:fit.nbasis(ne)
             temp = conv(thisPred,full(fit.bs{ne}(:,nb)));
             X = cat(2,X,temp(1:numel(y)));
         end
     end
     thisY = y(sum(X,2)~=0);
     thisX = X(sum(X,2)~=0,:);
     all_trials.firing = thisY;
     all_trials.design = thisX;
    
     [T(ni,:), y_pred] = czcrossvalstats(all_trials,model_parameters);
     if ni > 1
         fchange = T.fStat(ni)-T.fStat(ni-1);
         fprintf('Iter %d: fval = %.4f,  fchange = %.4f\n',ni,T.fStat(ni),fchange);
     else
         fprintf('Iter %d: fval = %.4f\n',ni,T.fStat(ni));
     end
     ni = ni+1;
end
fit.y_pred = y_pred; % last iteration's predicted data
fit.T = T;
fit.X = X;
fit.model_parameters = model_parameters; 
fit.kernelPreds = kernelPreds;
fit.gainPreds = gainPreds;
idx = 1;
for ne = 1:numel(fit.eventNames)
    fit.tstat{ne} = mdl2.Coefficients.tStat(idx:idx+numel(fit.gainNames{ne})+1-1);
    idx = idx+numel(fit.gainNames{ne})+1;
end