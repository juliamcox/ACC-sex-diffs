function plotRegression(plotVer, frameRate, savename,plotEvents,norm,histLen)

%plotVer: version of regression to plot (fHist to plot gcamp)
%frameRate
%savename: file name 
%norm: none, zscore, peak
%histLen: needs to be entered only if plotting fHist 


numPlots = 8;

fbasename = fullfile(whereAreWe('imaging'));

if contains(plotVer,'fHist')
    cons = 1:numel(plotEvents);
else
    [cons, con_shift, time_back_orig, time_forward_orig] = getEvents(plotVer,frameRate);
end

load(fullfile(fbasename,[savename '.mat']))

%% Plot heatmap of kernels for all neurons sorted by time of peak
[~,maxI_orig] = max(cell2mat(b_bs_all),[],2);
[~, maxI] = sort(maxI_orig);

%normalize regression coefficients 
if ~contains(plotVer,'fHist')
    switch norm
        case 'none'
            temp2 = cell2mat(b_bs_all);
            maxC = max(max(cell2mat(b_bs_all)));
            minC = min(min(cell2mat(b_bs_all)));
            minC = prctile(max(cell2mat(b_bs_all)),5);
            maxC = prctile(max(cell2mat(b_bs_all)),95); 
            if abs(minC)>abs(maxC)
                maxC = -minC;
            else
                minC = -maxC;
            end
        case 'zscore'
            temp2 = cell2mat(b_bs_all);
            mu = nanmean(temp2,2);
            sigma = nanstd(temp2,[],2);
            temp2 = (temp2-repmat(mu,1,size(temp2,2)))./repmat(sigma,1,size(temp2,2));
            maxC = 6;
            minC = -6;
        case 'peak'
            temp2 = cell2mat(b_bs_all);
            temp2=temp2./repmat(max((temp2),[],2),1,size(temp,2));
            maxC = 1;
            minC = -1;
    end
else
    maxC = max(max(max(cell2mat(b_bs_all))));
    minC = min(min(min(cell2mat(b_bs_all))));
    if abs(minC)>abs(maxC)
        maxC = -minC;
    else
        minC = -maxC;
    end
    temp2 = cell2mat(b_bs_all); 
end


screensize = get( groot, 'Screensize' );
screensize(4) = screensize(4)/1.5;
screensize(3) = screensize(3)/1.5;
figure('Position', screensize)

if ~contains(plotVer,'fHist')
    x=1;
    for ne = plotEvents
        idx = x:size(b_bs_all{ne},2)+x-1;
        a=subplot(1,numPlots,ne);
        if numel(time_forward_orig) > 1
            xaxis = linspace(-time_back_orig(ne),time_forward_orig(ne), size(b_bs_all{ne},2));
        else
            if con_shift(ne) == 0
                xaxis = linspace(-2,6,size(b_bs_all{ne},2));
            else
                xaxis = linspace(0,8,size(b_bs_all{ne},2));
            end
        end
        imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(maxI,idx), [minC maxC]);
        colormap(flipud(red2blue));

        box off
       % title(cons{ne});
        if ne == plotEvents(end)
            b=smallcolorbar(gca);
            b.Ticks = [minC 0 maxC];
            b.YLabel.String = 'Kernel value';
        end
        if ne ~= plotEvents(1)
            a.YTick = [];
        else
            ylabel('Neuron number (sorted by time of peak)')
        end
        x = x+size(b_bs_all{ne},2);
    end
    
else
    x=1;
    for ne = 1:numel(plotEvents)
        idx = x:size(b_bs_all{ne},2)+x-1;
        a=subplot(1,numPlots,ne);
        xaxis = linspace(histLen(ne,1),histLen(ne,2),size(b_bs_all{ne},2)); 
       
        imagesc(xaxis,1:size(b_bs_all{ne},1),temp2(maxI,idx),[minC maxC]);
        colormap(flipud(red2blue));
        box off
        title(plotEvents{ne});
        if ne == numel(plotEvents)
            b=smallcolorbar(gca);
           b.Ticks = [minC 0 maxC];
            b.YLabel.String = 'Kernel value';
        end
        if ne ~= 1
            a.YTick = [];
        else
            ylabel('Neuron number (sorted by time of peak)')
        end
        x = x+size(b_bs_all{ne},2);
    end
end
hold on, plot(xaxis,ones(size(xaxis)).*find(maxI_orig(maxI)>size(temp2,2)/numel(plotEvents),1,'first'),'k') 
