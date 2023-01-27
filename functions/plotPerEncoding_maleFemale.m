function [pval, chi2stat, posthoc] = plotPerEncoding_maleFemale(recs_f,recs_m,ver,frameRate, rawFlag, zscoreFlag,whichEvent)


plotParams = load(fullfile(whereAreWe('data'),'plotParams.mat'));
femaleC = plotParams.femaleC;
maleC = plotParams.maleC;
a = .01; %significance level

fname =  sprintf('linReg_fs%d_raw%d_zscore%d_basisSet_fig2_%s',frameRate,rawFlag,zscoreFlag,ver);
fbasename = fullfile(whereAreWe('imaging'));
fbasename_bs = fullfile(whereAreWe('basis_sets'));



[cons, ~, ~, ~, ~, ~, bsIDs] = getEvents(ver,frameRate); %get event names
b_female = [];
pmat_female = [];
for nr = 1:numel(recs_f)
    load(fullfile(fbasename,recs_f{nr},fname),'pvals','b','con_iden');
    load(fullfile(fbasename,recs_f{nr},'activeNeurons.mat'))
    b = cell2mat(b);
    b_female = cat(2, b_female,b(2:end,activeIdx));
    pvals = cell2mat(pvals');
    pmat = pvals(activeIdx,:)<a;
    pmat_female = cat(1,pmat_female,pmat);
end


pmat_male = [];
b_male = []; 
for nr = 1:numel(recs_m)
    load(fullfile(fbasename,recs_m{nr},fname),'pvals','b','con_iden');
    load(fullfile(fbasename,recs_m{nr},'activeNeurons.mat'))
    b = cell2mat(b);
    try
        b_male = cat(2, b_male,b(2:end,activeIdx));
    catch
        keyboard
    end
    pvals = cell2mat(pvals');
    pmat = pvals(activeIdx,:)<a;
    pmat_male = cat(1,pmat_male,pmat);
end

% Find kernels for each predictor and determine whether predominantly positive or negative
% load basis set

load(fullfile(fbasename_bs, ['bs_' bsIDs '.mat']))
bs = (full(eval(['bs_' bsIDs])));


% calculate kernel: extract coefficients for each event and multiply with basis set and sum
b_bs_all = cell(1,numel(cons));
for nn = 1:size(b_female,2)
    for ne = 1:numel(cons)
        thisWeights = b_female(con_iden==ne,nn);
        if iscell(bs)
            tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*bs{ne},2)';
        else
            tempWeights = sum(repmat(thisWeights',size(bs,1),1).*bs,2)';
        end
        b_bs_all{ne} = cat(1,b_bs_all{ne}, tempWeights);
        clear tempWeights
    end
end
b_bs_all_m = cell(1,numel(cons));
for nn = 1:size(b_male,2)
    for ne = 1:numel(cons)
        thisWeights = b_male(con_iden==ne,nn);
        if iscell(bs)
            tempWeights = sum(repmat(thisWeights',size(bs{ne},1),1).*bs{ne},2)';
        else
            tempWeights = sum(repmat(thisWeights',size(bs,1),1).*bs,2)';
        end
        b_bs_all_m{ne} = cat(1,b_bs_all_m{ne}, tempWeights);
        clear tempWeights
    end
end


% Find area under the curve for each kernel
for ne = whichEvent
    aoc_f(:,ne) = trapz(b_bs_all{ne}');
    aoc_m(:,ne) = trapz(b_bs_all_m{ne}');
    minD = min([aoc_f(:,ne); aoc_m(:,ne)]);
    maxD = max([aoc_f(:,ne); aoc_m(:,ne)]);
    figure()
    [h,b]=histcounts(aoc_f(:,ne),linspace(minD,maxD,30),'Normalization','probability');
    plot(b(1:end-1),h,'Color',femaleC)
    hold on
    [h,b]=histcounts(aoc_m(:,ne),linspace(minD,maxD,30),'Normalization','probability');
    plot(b(1:end-1),h,'Color',maleC)
    hold on, plot(mean(aoc_f(:,ne)), .2, 'v','Color',femaleC)
    hold on, plot(mean(aoc_m(:,ne)), .2, 'v','Color',maleC)
    title(cons{ne})
    box off
    set(gca,'FontSize',12)
    xlabel('Area under the curve')
    ylabel('Proportion of neurons')
end



for ne = whichEvent
    pos_f(ne) = sum(pmat_female(:,ne) & aoc_f(:,ne)>0);
    pos_m(ne) = sum(pmat_male(:,ne) & aoc_m(:,ne) > 0); 
    neg_f(ne) = sum(pmat_female(:,ne) & aoc_f(:,ne)<0);
    neg_m(ne) = sum(pmat_male(:,ne) & aoc_m(:,ne) < 0); 
    
    % Chi-squared test with post-hoc comparisons
    
    propTable(1,:) = [pos_f(ne) neg_f(ne) size(b_female,2)-pos_f(ne)-neg_f(ne)];
    propTable(2,:) = [pos_m(ne) neg_m(ne) size(b_male,2)-pos_m(ne)-neg_m(ne)];
    
    [pval(ne), chi2stat(ne), posthoc{ne}] = chi2PropTest_multi(propTable,a); 
    % relevant comparisons:
    % pos v neg female: 1v3, idx = 2
    % pos v neg male: 2v4, idx = 7
    % pos male v female: 1v2, idx = 1
    % neg male v female: 3v4, idx = 10
    
end

screensize = get( groot, 'Screensize' );
screensize(4) = screensize(4)/3.5;
screensize(3) = screensize(3)/1.5;
screensize(2) = screensize(2)+screensize(4)-100;
figure('Position', screensize)
for ne = whichEvent
    subplot(1,numel(cons),ne)
    b= bar([neg_f(ne)/size(b_female,2) neg_m(ne)/size(b_male,2);  pos_f(ne)/size(b_female,2) pos_m(ne)/size(b_male,2)],'EdgeColor','none');
    b(1).FaceColor = femaleC;
    b(2).FaceColor = maleC;
    box off
    ylabel('Proportion significant')
    set(gca,'XTickLabel',{'Pos';'Neg'})
    title(cons{ne})
end
