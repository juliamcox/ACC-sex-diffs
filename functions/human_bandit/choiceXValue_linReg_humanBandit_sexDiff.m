function stats=choiceXValue_linReg_humanBandit_sexDiff(choice_f,choice_m,bins,valType,zscoreFlag)



for nv = 1:numel(valType)
    
    nbins = numel(eval(sprintf('bins.%s;',valType{nv})));
    
    % extract data
    choice = [];
    sex = [];
    value = [];
    sIDs = [];
    trials = [];
    age = [];
    valueQuant = [];
    for ns = 1:size(choice_f.probRight_QDiff,1)
        ntrials = numel(eval(sprintf('choice_f.%s_trials{ns}',valType{nv})));
        thisvalue = eval(sprintf('choice_f.%s_trials{ns}',valType{nv}));
        thisbins = prctile(thisvalue,linspace(0,100,nbins+1));
        [~,~,thisvalue_quant] = histcounts(thisvalue,thisbins);
        choice = cat(1,choice,choice_f.choice_trials{ns});
        sex = cat(1,sex,ones(ntrials,1));
        value = cat(1,value,eval(sprintf('choice_f.%s_trials{ns}',valType{nv})));
        valueQuant = cat(1,valueQuant,thisvalue_quant);
        sIDs = cat(1,sIDs,ones(ntrials,1).*ns);
        age = cat(1,age,ones(ntrials,1).*choice_f.age(ns));
        trials = cat(1,trials,[1:ntrials]');
    end
    
    subCount = ns;
    
    for ns = 1:size(choice_m.probRight_QDiff,1)
        ntrials = numel(eval(sprintf('choice_m.%s_trials{ns}',valType{nv})));
        thisvalue = eval(sprintf('choice_m.%s_trials{ns}',valType{nv}));
        thisbins = prctile(thisvalue,linspace(0,100,nbins+1));
        [~,~,thisvalue_quant] = histcounts(thisvalue,thisbins);
        choice = cat(1,choice,choice_m.choice_trials{ns});
        sex = cat(1,sex,zeros(ntrials,1));
        value = cat(1,value,eval(sprintf('choice_m.%s_trials{ns}',valType{nv})));
        valueQuant = cat(1,valueQuant,thisvalue_quant);
        sIDs = cat(1,sIDs,ones(ntrials,1).*(ns+subCount));
        age = cat(1,age,ones(ntrials,1).*choice_m.age(ns));
        trials = cat(1,trials,[1:ntrials]');
    end
    
    
    
    % Value quantile w/o random slope
    X = table((choice),categorical(sex),categorical(valueQuant),categorical(sIDs),nanzscore(age),zscore(trials),'VariableNames',{'choice','sex','value','sID','age','trial'});
    
    f = 'choice ~ sex*value*age + (1+value|sID)';
    mdl = fitglme(X,f,'Distribution','binomial');
    eval(sprintf('stats.%s_quant = mdl;',valType{nv}));
    eval(sprintf('stats.%s_quant_coeff = dataset2cell(mdl.Coefficients);',valType{nv}));
    eval(sprintf('stats.%s_quant_anova = dataset2cell(anova(mdl));',valType{nv}));
    
    
    
    
end
