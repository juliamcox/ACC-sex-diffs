function fitLatencyDistributions_inverseGaussian(cohort,ext,zscoreFlag,savehere,laserType,valType_plot,groupingVar,ver,binNum,behaviorTable,latencyType)

% Fit 3-parameter shifted inverse gaussian distribution to response time data separated by group (sex, opsin) and groupingVar ('none','rew','value')
% ver: Parameterization 1 (matlab) or 2 


% Find opto mice 
ids_m = generateAnimalList(sprintf('%s_male',cohort));
ids_f = generateAnimalList(sprintf('%s_female',cohort));
ids_m_yfp = generateAnimalList(sprintf('%s_yfp_male',cohort));
ids_f_yfp = generateAnimalList(sprintf('%s_yfp_female',cohort));

idx = cell2mat(cellfun(@(x) find(strcmp(behaviorTable.aID,x)),cat(1,ids_m,ids_f,ids_m_yfp,ids_f_yfp),'UniformOutput',false));
behaviorTable = behaviorTable(idx,:); 
behaviorTable(behaviorTable.trial==1,:) = []; % remove first trial 

% Load plotting parameters
plotParams = load(fullfile(whereAreWe('figureCode'),'general_code','plotParams.mat'));

% Fitting parameters
thisValue = eval(sprintf('behaviorTable.%s_quant',valType_plot)); 
Latency   = eval(sprintf('behaviorTable.%s',latencyType)); 
Latency(Latency<0.01) = 0.01; % correct for different temporal resolutions 

behaviorTable.Latency = Latency;
behaviorTable.thisValue = thisValue;
behaviorTable(isnan(Latency),:) = [];

highCut = inf;
options.MaxIter = 2000; %increase number of iterations
options.MaxFunEvals = 2000; 
if ver == 1
    fits.thisParams = {'mu';'theta';'lambda'};
elseif ver == 2
    fits.thisParams = {'alpha';'theta';'gamma'};
end


% Which percentiles to plot for value 
minVal = min(thisValue)+1;
maxVal = max(thisValue)-1;


%% Fit custom inverse gaussian pdf

if ver == 1
    custpdf = @(data,mu,theta,lamda) pdf('InverseGaussian',data-theta,mu,lamda);
    custcdf = @(data,mu,theta,lamda) cdf('InverseGaussian',data-theta,mu,lamda);
elseif ver == 2
    custpdf = @(data,alpha,shift,gamma) pdf('InverseGaussian',data-shift,alpha./gamma,alpha.^2);
    custcdf = @(data,alpha,shift,gamma) cdf('InverseGaussian',data-shift,alpha./gamma,alpha.^2);
end

groups = {'f';'m';'f_yfp';'m_yfp'};

for ng = 1:numel(groups)
    clear thisFits
    thisIds = eval(sprintf('ids_%s',groups{ng}));
    for na = 1:numel(thisIds)
        switch groupingVar
            case 'none'
                if sum(behaviorTable.laser(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds{na}))==laserType) > 1
                    data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds{na})&(behaviorTable.laser==0)&behaviorTable.Latency<highCut);
                    try
                        thisFits{na}{1} = mle(data','pdf',custpdf,'start',[2 -.1 1],'lowerbound',[0 -1 0],'options',options);
                    catch
                        thisFits{na}{1} = mle(data','pdf',custpdf,'start',[.5 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                    end
                    
                    data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds{na})&(behaviorTable.laser==laserType)&behaviorTable.Latency<highCut);
                    try
                        thisFits{na}{2} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 1],'lowerbound',[0 -1 0],'options',options);
                    catch
                        thisFits{na}{2} = mle(data','pdf',custpdf,'start',[3 -.1 1],'lowerbound',[0 -1 0],'options',options);
                    end
                else
                    thisFits{na}{1} = [NaN NaN NaN];
                    thisFits{na}{2} = [NaN NaN NaN];
                end
            case 'rew'
                % Previously rewarded
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==0)&behaviorTable.previousReward==1);
                thisFits{1,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==laserType)&behaviorTable.previousReward==1);
                thisFits{1,na}{2} = mle(data','pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                % Previously unrewarded
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==0)&behaviorTable.previousReward==0);
                try
                    thisFits{2,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[2 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    thisFits{2,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                end
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==laserType)&behaviorTable.previousReward==0);
                
                try
                    thisFits{2,na}{2} = mle(data','pdf',custpdf,'cdf',custcdf,'start',[2 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    try
                        thisFits{2,na}{2} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[3 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                    catch
                        keyboard
                    end
                end
            case 'value'
                % High value
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==0)&behaviorTable.thisValue>=maxVal);
                try
                    thisFits{1,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    try
                    thisFits{1,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[2 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                    catch
                        keyboard
                    end
                end
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==laserType)&behaviorTable.thisValue>=maxVal);
                try
                    thisFits{1,na}{2} = mle(data','pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    thisFits{1,na}{2} = mle(data','pdf',custpdf,'start',[2 -.1 1],'lowerbound',[0 -1 0],'options',options);
                end
                % Low value
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==0)&behaviorTable.thisValue<=minVal);
                try
                    thisFits{2,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    thisFits{2,na}{1} = mle(data','cdf',custcdf,'pdf',custpdf,'start',[2 -.1 1],'lowerbound',[0 -1 0],'options',options);
                end
                data = behaviorTable.Latency(behaviorTable.laserSession==1&strcmp(behaviorTable.aID,thisIds(na))&(behaviorTable.laser==laserType)&behaviorTable.thisValue<=minVal);
                try
                    thisFits{2,na}{2} = mle(data','pdf',custpdf,'start',[1 -.1 .5],'lowerbound',[0 -1 0],'options',options);
                catch
                    thisFits{2,na}{2} = mle(data','pdf',custpdf,'start',[2 -.1 1],'lowerbound',[0 -1 0],'options',options);
                end
                
        end
    end
    eval(sprintf('fits.fits_%s = thisFits;',groups{ng}))
    clear thisFits
end



fits.ids_f = ids_f;
fits.ids_m = ids_m;
fits.ids_f_yfp = ids_f_yfp;
fits.ids_m_yfp = ids_m_yfp;

save(fullfile(savehere,sprintf('InverseGaussian_%s_%s_zscore%d_grouping_%s_ver%d_laserType%d.mat',cohort,ext,zscoreFlag,groupingVar,ver,laserType)),'fits')

