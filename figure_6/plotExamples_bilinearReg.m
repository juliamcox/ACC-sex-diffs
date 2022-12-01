function plotExamples_bilinearReg(pvals,alpha,ver,thisFnames)

ptiles = linspace(0,100,4);
rawFlag = 1;
frameRate = 20; 
sigMat = pvals<alpha;
neurons = 1:size(sigMat,1);

% saveLoc = fullfile(whereAreWe('figurecode'),'processed_data','fig6',ver);
% if ~exist(saveLoc,'dir')
%     mkdir(saveLoc);
% end

fbasename = fullfile(whereAreWe('bilinear'),'dataset_single1');

figure('Position',[440 105 466 711]);
fh(1) = gcf; set(fh(1),'Visible','off')
figure('Position',[311 296 907 420]);
fh(2) = gcf; set(fh(2),'Visible','off')



for nr = 1:numel(neurons)
    load(fullfile(fbasename,ver,thisFnames{neurons(nr)}))
    if rawFlag == 0
        temp = thisFnames{neurons(nr)};
        idx = strfind(temp,'rawFlag0');
        temp(idx:idx+7) = 'rawFlag1'; 
        load(fullfile(fbasename,'rawFlag1_frameRate20_qLearn_session_all',temp));
    else
        load(fullfile(fbasename,'rawFlag1_frameRate20_qLearn_session_all',thisFnames{neurons(nr)}));
    end
    y = nanzscore(y);
    fit_real = extractEvents_bilinear(fit_real,y,ptiles);
    for n = 1:numel(fit_crossVal.T)
        r2(n) = fit_crossVal.T{n}.R2(end);
    end
    r2=mean(r2);
    savename = fullfile(sprintf('%s_neuron%d_f',ver,neurons(nr)));
    plotBilinearRegression(fit_real,r2,savename,saveLoc,fh);
    %close all
    clear fit_real y
end