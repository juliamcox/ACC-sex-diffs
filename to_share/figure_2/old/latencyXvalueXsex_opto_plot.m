function stats = latencyXvalueXsex_opto_plot(latency_f, latency_m, latency_f_yfp, latency_m_yfp, epochs,valType,latencyType,zscoreFlag,bins)

%%% Figure 4

% Plot latency x value x sex as extracted in latencyXvalueXsex_opto.m

% ext: file extension (LAS,UNI et)
% epochs: laser conditions to plot
% valType: qDiff, qTot, qChosen, qUnchosen, qChosenDiff
% latencyType: which latency to plot
% zscoreFlag: 0 or 1 to zscore resposnse times
% cohort: basename for cohort (e.g. 'ACC_DMS_nphr') for animal lists generated with generateAnimalList
% bins 

plotParams = load(fullfile(whereAreWe('figureCode'),'general_code','plotParams.mat')); 

markersize = 10;


%% Plot trial initiation latency for males and females as a function of value (quantiles)


for nv = 1:numel(valType)
    for nl = 1:numel(latencyType)
        for ne = 1:numel(epochs)
            % Plot means only
           
            
            groupvec = [];
            vec = [];
            for na = 1:size(mu_f,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_f(na,1:end-1)');
            end
            for na = 1:size(mu_m,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_m(na,1:end-1)');
            end
            for na = 1:size(mu_f_laser,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_f_laser(na,1:end-1)');
            end
            for na = 1:size(mu_m_laser,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_m_laser(na,1:end-1)');
            end
            
            
            
            for na = 1:size(mu_f_yfp,1)
                groupvec =cat(1,groupvec, [[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_f_yfp(na,1:end-1)');
            end
            for na = 1:size(mu_m_yfp,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))']);
                vec =cat(1,vec,mu_m_yfp(na,1:end-1)');
            end
            for na = 1:size(mu_f_laser_yfp,1)
                groupvec =cat(1,groupvec, [[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))']);
                vec = cat(1,vec,mu_f_laser_yfp(na,1:end-1)');
            end
            for na = 1:size(mu_m_laser_yfp,1)
                groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))']);
                vec =cat(1,vec,mu_m_laser_yfp(na,1:end-1)');
            end
            
            [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec,'varnames', {'value';'sex';'opsin';'laser'}, 'display', 'off', 'model','full');
            try
                tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3 4],'CType','hsd','Display','off');
            catch
                try
                    [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec(:,[1 3 4]),'varnames', {'value';'opsin';'laser'}, 'display', 'off', 'model','full');
                    tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3],'CType','hsd','Display','off');
                catch
                    [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec(:,[1 2 4]),'varnames', {'value';'sex';'laser'}, 'display', 'off', 'model','full');
                    tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3],'CType','hsd','Display','off');
                end
            end
            eval(sprintf('stats.%s_%s_%s_quant = tempStats;', valType{nv}, latencyType{nl}, epochs{ne}))
            
           %%
              
              pause(.5); f=figure('Units','inches','Position',[5,5,9 4]);
              thisbin = eval(sprintf('bins.%s', valType{nv}));
              
              % Plot female laser and no laser trials
              subplot(2,2,1); hold on
              xaxis = thisbin+mean(diff(thisbin))/2;
              
              mu_f = eval(sprintf('latency_f.%s_%s_quant;',latencyType{nl},valType{nv}));
              if isempty(mu_f)
                  mu_f = nan(1,length(xaxis));
              end
              mu_f_laser = eval(sprintf('latency_f.%s_%s_%s_quant',latencyType{nl},valType{nv},epochs{ne}));
              if isempty(mu_f_laser)
                  mu_f_laser = nan(1,length(xaxis));
              end
              mu_f = mu_f(~isnan(mu_f_laser(:,1)),:);
              mu_f_laser = mu_f_laser(~isnan(mu_f_laser(:,1)),:);
              sem_f_laser = nanstd(mu_f_laser,[],1)./sqrt(size(mu_f_laser,1));
              sem_f = nanstd(mu_f,[],1)./sqrt(size(mu_f,1));
              
              b=bar(xaxis,[nanmean(mu_f,1);nanmean(mu_f_laser,1)]');
              pause(.01);
              scatter(xaxis, mean(mu_f),15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor','none')
              scatter(xaxis, nanmean(mu_f_laser),15, 'MarkerFaceColor', plotParams.femaleC, 'MarkerEdgeColor','none')
              
              b(1).FaceColor = 'none';
              b(2).FaceColor = 'none';
              b(1).EdgeColor = 'none';
              b(2).EdgeColor = 'none';
              b(1).LineWidth = 1.5;
              b(2).LineWidth = 1.5;
              
              for nb = 1:numel(xaxis)
                  tempX = repmat(xaxis(nb)+b(1).XOffset,1,size(mu_f,1));%+datasample(-.05:.01:.05,size(mu_f,1));
                  scatter(tempX,mu_f(:,nb),15,'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              
              for nb = 1:numel(xaxis)
                  tempX = repmat(xaxis(nb)+b(2).XOffset,1,size(mu_f,1));%+datasample(-.05:.01:.05,size(mu_f,1));
                  scatter(tempX,mu_f_laser(:,nb),15,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              legendh=errorbar(xaxis,nanmean(mu_f,1), sem_f, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
              
              legendh(2)=errorbar(xaxis,nanmean(mu_f_laser,1), sem_f_laser, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.femaleC 1],'LineWidth',1.5);

%
%               plot(repmat(xaxis,size(mu_f,1),1)', mu_f','-','Color', [.3 .3 .3 .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_f,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.3 .3 .3],'MarkerFaceAlpha',1, 'LineWidth', .5);
%               
%               plot(repmat(xaxis,size(mu_f_laser,1),1)',mu_f_laser','-','Color',[plotParams.femaleC .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_f_laser,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotParams.femaleC,'MarkerFaceAlpha',1);
%                legendh=errorbar(xaxis,nanmean(mu_f,1), sem_f, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
%               
%               legendh(2)=errorbar(xaxis,nanmean(mu_f_laser,1), sem_f_laser, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.femaleC 1],'LineWidth',1.5);
%               
              axis square
              
              
              title('Female: NpHR')
              set(gca,'FontSize',12,'XLim', [thisbin(1)-.25 thisbin(end)+.25])
              switch latencyType{nl}
                  case 'trialStart'
                    if zscoreFlag == 1
                        ylabel('Trial initiation latency (zscore)')
                    elseif zscoreFlag == 2
                        ylabel('Trial initiation latency (log(sec))')
                    else
                        ylabel('Trial initiation latency (sec)')
                        set(gca,'YLim',[0 20],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                    end
                  case 'leverPress'
                      if zscoreFlag == 1
                          ylabel('Lever press latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Lever press latency (log(sec))')
                      else
                          ylabel('Lever press latency (sec)')
                          set(gca,'YLim',[0 3],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                      end
                  case 'withdraw'
                      if zscoreFlag == 1
                          ylabel('Nose poke withdrawal latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Nose poke withdrawal latency (log(sec))')
                      else
                          ylabel('Nose poke withdrawal latency (sec)')
                          set(gca,'YLim',[0 .25],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                          
                      end
              end
              xlabel(valType{nv})
              if zscoreFlag==1
                  set(gca,'YLim',[-.25 .36],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              elseif zscoreFlag == 2
                  set(gca,'YLim',[-2 1],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              end
              
              
              % Plot male laser and no laser trials
              subplot(2,2,2); hold on
              mu_m = eval(sprintf('latency_m.%s_%s_quant;',latencyType{nl},valType{nv}));
              if isempty(mu_m)
                  mu_m = nan(1,length(xaxis));
              end
              mu_m_laser = eval(sprintf('latency_m.%s_%s_%s_quant',latencyType{nl},valType{nv},epochs{ne}));
              if isempty(mu_m_laser)
                  mu_m_laser = nan(1,length(xaxis));
              end
              
              mu_m = mu_m(~isnan(mu_m_laser(:,1)),:);
              mu_m_laser = mu_m_laser(~isnan(mu_m_laser(:,1)),:);
              sem_m_laser = nanstd(mu_m_laser,[],1)./sqrt(size(mu_m_laser,1));
              sem_m = nanstd(mu_m)./sqrt(size(mu_m,1));
              
              b=bar(xaxis,[nanmean(mu_m,1);nanmean(mu_m_laser,1)]');
              
              pause(.01); 
               scatter(xaxis, mean(mu_m),15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor','none')
              scatter(xaxis, nanmean(mu_m_laser),15, 'MarkerFaceColor', plotParams.maleC, 'MarkerEdgeColor','none')
            
              b(1).FaceColor = 'none';
              b(2).FaceColor = 'none';
              b(1).EdgeColor = 'none';
              b(2).EdgeColor = 'none';
              b(1).LineWidth = 1.5;
              b(2).LineWidth = 1.5; 
              
              for nb = 1:numel(xaxis)
                 tempX = repmat(xaxis(nb)+b(1).XOffset,1,size(mu_m,1));%+datasample(-.05:.01:.05,size(mu_m,1));
                 scatter(tempX,mu_m(:,nb),15,'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              
              for nb = 1:numel(xaxis)
                 tempX = repmat(xaxis(nb)+b(2).XOffset,1,size(mu_m,1));%+datasample(-.05:.01:.05,size(mu_m,1));
                 scatter(tempX,mu_m_laser(:,nb),15,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              legendh=errorbar(xaxis,nanmean(mu_m,1), sem_m, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
              
              legendh(2)=errorbar(xaxis,nanmean(mu_m_laser,1), sem_m_laser, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.maleC 1],'LineWidth',1.5);
              
              
%               plot(repmat(xaxis,size(mu_m,1),1)', mu_m', '-','Color',[.3 .3 .3 .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_m,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.3 .3 .3],'MarkerFaceAlpha',1);
%               
%               plot(repmat(xaxis,size(mu_m_laser,1),1)', mu_m_laser','-','Color',[plotParams.maleC .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_m_laser,1), markersize,  'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotParams.maleC,'MarkerFaceAlpha',1);
%               legendh(1)=errorbar(xaxis,nanmean(mu_m,1), sem_m, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3],'LineWidth',1.5);
%               legendh(2)=errorbar(xaxis,nanmean(mu_m_laser,1), sem_m_laser, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.maleC 1],'LineWidth',1.5);
%               
              title('Male: NpHR')
              set(gca,'FontSize',12,'XLim', [thisbin(1)-.25 thisbin(end)+.25])
              switch latencyType{nl}
                  case 'trialStart'
                    if zscoreFlag == 1
                        ylabel('Trial initiation latency (zscore)')
                    elseif zscoreFlag == 2
                        ylabel('Trial initiation latency (log(sec))')
                    else
                        ylabel('Trial initiation latency (sec)')
                        set(gca,'YLim',[0 20],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                    end
                  case 'leverPress'
                      if zscoreFlag == 1
                          ylabel('Lever press latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Lever press latency (log(sec))')
                      else
                          ylabel('Lever press latency (sec)')
                          set(gca,'YLim',[0 3],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                      end
                  case 'withdraw'
                      if zscoreFlag == 1
                          ylabel('Nose poke withdrawal latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Nose poke withdrawal latency (log(sec))')
                      else
                          ylabel('Nose poke withdrawal latency (sec)')
                          set(gca,'YLim',[0 .25],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                          
                      end
              end
              xlabel(valType{nv})
              if zscoreFlag==1
                  set(gca,'YLim',[-.25 .36],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              elseif zscoreFlag == 2
                  set(gca,'YLim',[-2 1],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              end
              
              axis square
              
              % Plot YFP
              subplot(2,2,3); hold on
              
              mu_f_yfp = eval(sprintf('latency_f_yfp.%s_%s_quant;',latencyType{nl},valType{nv}));
              if isempty(mu_f_yfp)
                  mu_f_yfp = nan(1,length(xaxis));
              end
              mu_f_laser_yfp = eval(sprintf('latency_f_yfp.%s_%s_%s_quant',latencyType{nl},valType{nv},epochs{ne}));
              if isempty(mu_f_laser_yfp)
                  mu_f_laser_yfp = nan(1,length(xaxis));
              end
              
              mu_f_yfp = mu_f_yfp(~isnan(mu_f_laser_yfp(:,1)),:);
              mu_f_laser_yfp = mu_f_laser_yfp(~isnan(mu_f_laser_yfp(:,1)),:);
              sem_f_laser_yfp = nanstd(mu_f_laser_yfp,[],1)./sqrt(size(mu_f_laser_yfp,1));
              sem_f_yfp = nanstd(mu_f_yfp,[],1)./sqrt(size(mu_f_yfp,1));
              
              b=bar(xaxis,[nanmean(mu_f_yfp,1);nanmean(mu_f_laser_yfp,1)]');
                            pause(.01); 

              scatter(xaxis, mean(mu_f_yfp),15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor','none')
              scatter(xaxis, nanmean(mu_f_laser_yfp),15, 'MarkerFaceColor', plotParams.femaleC, 'MarkerEdgeColor','none')
            
              b(1).FaceColor = 'none';
              b(2).FaceColor = 'none';
              b(1).EdgeColor = 'none';
              b(2).EdgeColor = 'none';
              b(1).LineWidth = 1.5;
              b(2).LineWidth = 1.5; 
              
              for nb = 1:numel(xaxis)
                 tempX = repmat(xaxis(nb)+b(1).XOffset,1,size(mu_f_yfp,1));%+datasample(-.05:.01:.05,size(mu_f_yfp,1));
                 scatter(tempX,mu_f_yfp(:,nb),15,'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              
              for nb = 1:numel(xaxis)
                 tempX = repmat(xaxis(nb)+b(2).XOffset,1,size(mu_f_yfp,1));%+datasample(-.05:.01:.05,size(mu_f_yfp,1));
                 scatter(tempX,mu_f_laser_yfp(:,nb),15,'MarkerFaceColor',plotParams.femaleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              legendh=errorbar(xaxis,nanmean(mu_f_yfp,1), sem_f_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
              
              legendh(2)=errorbar(xaxis,nanmean(mu_f_laser_yfp,1), sem_f_laser_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.femaleC 1],'LineWidth',1.5);
              
              
              
              
%               plot(repmat(xaxis,size(mu_f_yfp,1),1)', mu_f_yfp','-','Color',[.3 .3 .3 .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_f_yfp,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.3 .3 .3],'MarkerFaceAlpha',1);
%               
%               plot(repmat(xaxis,size(mu_f_laser_yfp,1),1)', mu_f_laser_yfp','-','Color',[plotParams.femaleC .6], 'LineWidth', .5);
%               scatter(xaxis,nanmean(mu_f_laser_yfp,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotParams.femaleC,'MarkerFaceAlpha',1);
%               errorbar(xaxis,nanmean(mu_f_yfp,1), sem_f_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3],'LineWidth',1.5);
%               errorbar(xaxis,nanmean(mu_f_laser_yfp,1), sem_f_laser_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.femaleC 1],'LineWidth',1.5);
%               
              title('Female: YFP')
              set(gca,'FontSize',12,'XLim', [thisbin(1)-.25 thisbin(end)+.25])
              switch latencyType{nl}
                  case 'trialStart'
                    if zscoreFlag == 1
                        ylabel('Trial initiation latency (zscore)')
                    elseif zscoreFlag == 2
                        ylabel('Trial initiation latency (log(sec))')
                    else
                        ylabel('Trial initiation latency (sec)')
                        set(gca,'YLim',[0 20],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                    end
                  case 'leverPress'
                      if zscoreFlag == 1
                          ylabel('Lever press latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Lever press latency (log(sec))')
                      else
                          ylabel('Lever press latency (sec)')
                          set(gca,'YLim',[0 3],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                      end
                  case 'withdraw'
                      if zscoreFlag == 1
                          ylabel('Nose poke withdrawal latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Nose poke withdrawal latency (log(sec))')
                      else
                          ylabel('Nose poke withdrawal latency (sec)')
                          set(gca,'YLim',[0 .25],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                          
                      end
              end
              xlabel(valType{nv})
              if zscoreFlag==1
                  set(gca,'YLim',[-.25 .36],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              elseif zscoreFlag == 2
                  set(gca,'YLim',[-2 1],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              end
              
              axis square
              
              subplot(2,2,4); hold on
              mu_m_yfp = eval(sprintf('latency_m_yfp.%s_%s_quant;',latencyType{nl},valType{nv}));
              if isempty(mu_m_yfp)
                  mu_m_yfp = nan(1,length(xaxis));
              end
              mu_m_laser_yfp = eval(sprintf('latency_m_yfp.%s_%s_%s_quant',latencyType{nl},valType{nv},epochs{ne}));
              if isempty(mu_m_laser_yfp)
                  mu_m_laser_yfp = nan(1,length(xaxis));
              end
              
              mu_m_yfp = mu_m_yfp(~isnan(mu_m_laser_yfp(:,1)),:);
              mu_m_laser_yfp = mu_m_laser_yfp(~isnan(mu_m_laser_yfp(:,1)),:);
              sem_m_laser_yfp = nanstd(mu_m_laser_yfp,[],1)./sqrt(sum(~isnan(mu_m_laser_yfp(:,1))));
              sem_m_yfp = nanstd(mu_m_yfp,[],1)./sqrt(size(mu_m_yfp,1));

              
              %               plot(repmat(xaxis,size(mu_m_yfp,1),1)', mu_m_yfp','-','Color',[.3 .3 .3 .6], 'LineWidth', .5);
              %               scatter(xaxis,nanmean(mu_m_yfp,1), markersize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.3 .3 .3],'MarkerFaceAlpha',1);
              %
              %               plot(repmat(xaxis,size(mu_m_laser_yfp,1),1)', mu_m_laser_yfp','-','Color',[plotParams.maleC .6], 'LineWidth', .5);
              %               scatter(xaxis,nanmean(mu_m_laser_yfp,1),markersize,  'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotParams.maleC,'MarkerFaceAlpha',1);
              %               legendh=errorbar(xaxis,nanmean(mu_m_yfp,1), sem_m_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3],'LineWidth',1.5);
              %               legendh(2)=errorbar(xaxis,nanmean(mu_m_laser_yfp,1), sem_m_laser_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.maleC 1],'LineWidth',1.5);
              %
              
              b=bar(xaxis,[nanmean(mu_m_yfp,1);nanmean(mu_m_laser_yfp,1)]');
              
              pause(.001);
              scatter(xaxis, mean(mu_m_yfp),15, 'MarkerFaceColor', [.3 .3 .3], 'MarkerEdgeColor','none')
              scatter(xaxis, nanmean(mu_m_laser_yfp),15, 'MarkerFaceColor', plotParams.maleC, 'MarkerEdgeColor','none')
              
              b(1).FaceColor = 'none';
              b(2).FaceColor = 'none';
              b(1).EdgeColor = 'none';
              b(2).EdgeColor = 'none';
              b(1).LineWidth = 1.5;
              b(2).LineWidth = 1.5;
              
              for nb = 1:numel(xaxis)
                  tempX = repmat(xaxis(nb)+b(1).XOffset,1,size(mu_m_yfp,1));%+datasample(-.05:.01:.05,size(mu_m_yfp,1));
                  scatter(tempX,mu_m_yfp(:,nb),15,'MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              
              for nb = 1:numel(xaxis)
                  tempX = repmat(xaxis(nb)+b(2).XOffset,1,size(mu_m_yfp,1));%+datasample(-.05:.01:.05,size(mu_m_yfp,1));
                  scatter(tempX,mu_m_laser_yfp(:,nb),15,'MarkerFaceColor',plotParams.maleC,'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
              end
              legendh=errorbar(xaxis,nanmean(mu_m_yfp,1), sem_m_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [.3 .3 .3], 'LineWidth', 1.5);
              
              legendh(2)=errorbar(xaxis,nanmean(mu_m_laser_yfp,1), sem_m_laser_yfp, 'LineStyle', '-', 'CapSize', 0, 'Color', [plotParams.maleC 1],'LineWidth',1.5);

              title('Male: YFP')
              set(gca,'FontSize',12,'XLim', [thisbin(1)-.25 thisbin(end)+.25])
              switch latencyType{nl}
                  case 'trialStart'
                    if zscoreFlag == 1
                        ylabel('Trial initiation latency (zscore)')
                    elseif zscoreFlag == 2
                        ylabel('Trial initiation latency (log(sec))')
                    else
                        ylabel('Trial initiation latency (sec)')
                        set(gca,'YLim',[0 20],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                    end
                  case 'leverPress'
                      if zscoreFlag == 1
                          ylabel('Lever press latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Lever press latency (log(sec))')
                      else
                          ylabel('Lever press latency (sec)')
                          set(gca,'YLim',[0 3],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                      end
                  case 'withdraw'
                      if zscoreFlag == 1
                          ylabel('Nose poke withdrawal latency (zscore)')
                      elseif zscoreFlag == 2
                          ylabel('Nose poke withdrawal latency (log(sec))')
                      else
                          ylabel('Nose poke withdrawal latency (sec)')
                          set(gca,'YLim',[0 .25],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
                          
                      end
              end
              xlabel(valType{nv})
              if zscoreFlag==1
                  set(gca,'YLim',[-.25 .36],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              elseif zscoreFlag == 2
                  set(gca,'YLim',[-2 1],'XTick',linspace(thisbin(1),thisbin(end),length(xaxis)-1))
              end
              
              axis square
              
              
              acounter = 1;
              groupvec = [];
              vec = [];
              for na = 1:size(mu_f,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_f(na,1:end-1)');
                  acounter = acounter+1;
              end
              for na = 1:size(mu_m,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_m(na,1:end-1)');
                  acounter = acounter+1;
              end
              acounter=1;
              for na = 1:size(mu_f_laser,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_f_laser(na,1:end-1)');
                  acounter=acounter+1;
              end
              for na = 1:size(mu_m_laser,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_m_laser(na,1:end-1)');
                  acounter=acounter+1;
              end
              
              
              acounter_init = acounter;
              for na = 1:size(mu_f_yfp,1)
                  groupvec =cat(1,groupvec, [[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_f_yfp(na,1:end-1)');
                  acounter = acounter+1;
              end
              for na = 1:size(mu_m_yfp,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec =cat(1,vec,mu_m_yfp(na,1:end-1)');
                  acounter=acounter+1; 
              end
              acounter=acounter_init;
              for na = 1:size(mu_f_laser_yfp,1)
                  groupvec =cat(1,groupvec, [[1:numel(thisbin(1:end-1))]'  zeros(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec = cat(1,vec,mu_f_laser_yfp(na,1:end-1)');
                  acounter=acounter+1;
              end
              for na = 1:size(mu_m_laser_yfp,1)
                  groupvec = cat(1,groupvec,[[1:numel(thisbin(1:end-1))]'  ones(size(thisbin(1:end-1)))' zeros(size(thisbin(1:end-1)))' ones(size(thisbin(1:end-1)))' acounter.*ones(size(thisbin(1:end-1)))']);
                  vec =cat(1,vec,mu_m_laser_yfp(na,1:end-1)');
                  acounter=acounter+1;
              end
              
              [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec(:,1:4),'varnames', {'value';'sex';'opsin';'laser'}, 'display', 'off', 'model','full');
              try
                  tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3 4],'CType','hsd','Display','off');
              catch
                  try
                      [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec(:,[1 3 4]),'varnames', {'value';'opsin';'laser'}, 'display', 'off', 'model','full');
                      tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3],'CType','hsd','Display','off');
                  catch
                      [tempStats.p_anova,tempStats.t_anova,tempStats.s_anova] = anovan(vec, groupvec(:,[1 2 4]),'varnames', {'value';'sex';'laser'}, 'display', 'off', 'model','full');
                      tempStats.postHoc_tukey = multcompare(tempStats.s_anova,'dimension',[1 2 3],'CType','hsd','Display','off');
                  end
              end
              eval(sprintf('stats.%s_%s_%s_quant = tempStats;', valType{nv}, latencyType{nl}, epochs{ne}))
              
              
              % Post hoc tests: laser v no laser for each value per group
              % test for normality 
           
              for nb = 1:size(mu_f,2)-1
                  if lillietest(mu_f(:,nb)) || lillietest(mu_f_laser(:,nb))
                      [stats.quant_posthoc_female_laserNoLaser(nb),~] = signrank(mu_f(:,nb),mu_f_laser(:,nb));
                      stats.quant_posthoc_female_laserNoLaser_ttest(nb) = 0;
                  else
                      [~,stats.quant_posthoc_female_laserNoLaser(nb)] = ttest(mu_f(:,nb),mu_f_laser(:,nb));
                      stats.quant_posthoc_female_laserNoLaser_ttest(nb) = 1;
                  end
              end
              
              for nb = 1:size(mu_m,2)-1
                  if lillietest(mu_m(:,nb)) || lillietest(mu_m_laser(:,nb))
                      [stats.quant_posthoc_male_laserNoLaser(nb),~] = signrank(mu_m(:,nb),mu_m_laser(:,nb));
                      stats.quant_posthoc_male_laserNoLaser_ttest(nb) = 0;
                  else
                      [~,stats.quant_posthoc_male_laserNoLaser(nb)] = ttest(mu_m(:,nb),mu_m_laser(:,nb));
                      stats.quant_posthoc_male_laserNoLaser_ttest(nb) = 1;
                  end
              end
              
              for nb = 1:size(mu_f_yfp,2)-1
                  if lillietest(mu_f_yfp(:,nb)) || lillietest(mu_f_laser(:,nb))
                      [stats.quant_posthoc_female_laserNoLaser_yfp(nb),~] = signrank(mu_f_yfp(:,nb),mu_f_laser_yfp(:,nb));
                      stats.quant_posthoc_female_laserNoLaser_yfp_ttest(nb) = 0;
                  else
                      [~,stats.quant_posthoc_female_laserNoLaser_yfp(nb)] = ttest(mu_f_yfp(:,nb),mu_f_laser_yfp(:,nb));
                      stats.quant_posthoc_female_laserNoLaser_yfp_ttest(nb) = 1;
                  end
              end
              
              for nb = 1:size(mu_m_yfp,2)-1
                  if lillietest(mu_m_yfp(:,nb)) || lillietest(mu_m_laser_yfp(:,nb))
                      [stats.quant_posthoc_male_laserNoLaser_yfp(nb),~] = signrank(mu_m_yfp(:,nb),mu_m_laser_yfp(:,nb));
                      stats.quant_posthoc_male_laserNoLaser_yfp_ttest(nb) = 0;
                  else
                      [~,stats.quant_posthoc_male_laserNoLaser_yfp(nb)] = ttest(mu_m_yfp(:,nb),mu_m_laser_yfp(:,nb));
                      stats.quant_posthoc_male_laserNoLaser_yfp_ttest(nb) = 1;
                  end
              end
              
             % YFP v Opsin for laser 
             for nb = 1:size(mu_f,2)-1
                 if lillietest(mu_f_laser_yfp(:,nb)) || lillietest(mu_f_laser(:,nb))
                     [stats.quant_posthoc_female_laser_opsin(nb),~] = ranksum(mu_f_laser(:,nb),mu_f_laser_yfp(:,nb));
                     stats.quant_posthoc_female_laser_opsin_ttest(nb) = 0;
                 else
                     [~,stats.quant_posthoc_female_laser_opsin(nb)]= ttest2(mu_f_laser(:,nb),mu_f_laser_yfp(:,nb));
                     stats.quant_posthoc_female_laser_opsin_ttest(nb) = 1;
                 end
             end
             
             for nb = 1:size(mu_m,2)-1
                 if lillietest(mu_m_laser_yfp(:,nb)) || lillietest(mu_m_laser(:,nb))
                     [stats.quant_posthoc_male_laser_opsin(nb),~] = ranksum(mu_m_laser(:,nb),mu_m_laser_yfp(:,nb));
                     stats.quant_posthoc_male_laser_opsin_ttest(nb) = 0;
                 else
                     [~,stats.quant_posthoc_male_laser_opsin(nb)]= ttest2(mu_m_laser(:,nb),mu_m_laser_yfp(:,nb));
                     stats.quant_posthoc_male_laser_opsin_ttest(nb) = 1;
                 end
             end
              
              
              
              tbl = table(vec);
              sex = nominal(groupvec(:,2));
              value = nominal(groupvec(:,1));
              tbl.sex = sex;
              tbl.value = value;
              subject = nominal(groupvec(:,5));
              tbl.subject = subject;
              laser = nominal(groupvec(:,4));
              opsin = nominal(groupvec(:,3)); 
              tbl.laser = laser;
              tbl.opsin = opsin;
              stats.mdl_quant = fitlme(tbl, 'vec~value*sex*opsin*laser+(1|subject)','DummyVarCoding','effects');
          

          end
      end
  end

