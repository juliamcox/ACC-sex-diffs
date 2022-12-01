function stats=latencyXValue_linReg_humanBandit_sexDiff(latency_f,latency_m,bins,valType,latencyType,zscoreFlag)

load(fullfile(whereAreWe('bucket'), 'Manuscript_figures','plotParams.mat'))
% Cameron's thesis colors 
% femaleC = [248, 118, 109]./255;
% maleC = [0, 191, 196]./255;

for nl = 1:numel(latencyType)
    for nv = 1:numel(valType)
        
        thisbins = eval(sprintf('bins.%s;',valType{nv})); 
        thisbins = thisbins+mean(diff(thisbins))/2; % center bins
   
        % extract data
        latency = [];
        sex = [];
        value = []; 
        sIDs = [];
        trials = [];
        age = [];    
        valueQuant = [];
        for ns = 1:size(latency_f.choice_QDiff,1)
            ntrials = numel(eval(sprintf('latency_f.%s_trials{ns}',latencyType{nv})));
            latency = cat(1,latency,eval(sprintf('latency_f.%s_trials{ns}',latencyType{nv})));
            sex = cat(1,sex,ones(ntrials,1));
            value = cat(1,value,eval(sprintf('latency_f.%s_value{ns}',valType{nv})));
            valueQuant = cat(1,valueQuant,eval(sprintf('latency_f.%s_valueQuant{ns}',valType{nv})));
            sIDs = cat(1,sIDs,ones(ntrials,1).*ns);
            age = cat(1,age,ones(ntrials,1).*latency_f.age(ns));
            trials = cat(1,trials,[1:ntrials]');
        end
        
        subCount = ns; 
        
        for ns = 1:size(latency_m.choice_QDiff,1)
            ntrials = numel(eval(sprintf('latency_m.%s_trials{ns}',latencyType{nv})));
            latency = cat(1,latency,eval(sprintf('latency_m.%s_trials{ns}',latencyType{nv})));
            sex = cat(1,sex,zeros(ntrials,1));
            value = cat(1,value,eval(sprintf('latency_m.%s_value{ns}',valType{nv})));
            valueQuant = cat(1,valueQuant,eval(sprintf('latency_m.%s_valueQuant{ns}',valType{nv})));
            sIDs = cat(1,sIDs,ones(ntrials,1).*(ns+subCount)); 
            age = cat(1,age,ones(ntrials,1).*latency_m.age(ns));
            trials = cat(1,trials,[1:ntrials]');
        end
        
        latency(latency==0) = NaN;
        X = table((latency),categorical(sex),zscore(value),categorical(sIDs),nanzscore(age),zscore(trials),categorical(valueQuant),'VariableNames',{'latency','sex','value','sID','age','trial','valueQuant'});
        

       
        
        f = 'latency ~ sex*valueQuant*age+ (1|sID)';
        mdl = fitlme(X,f,'DummyVarCoding','effect');
        eval(sprintf('stats.%s_quant_quant = mdl;',valType{nv}));
        eval(sprintf('stats.%s_quant_coeff = dataset2cell(mdl.Coefficients);',valType{nv}));
        eval(sprintf('stats.%s_quant_anova = dataset2cell(anova(mdl,''DFMethod'',''Satterthwaite''));',valType{nv}));
 

%         eval(sprintf('stats.%s_quant_contrasts(1) = coefTest(mdl, [0 -2 0 0 0 0 0 -2 0 0 0]);', valType{nv}));
%         eval(sprintf('stats.%s_quant_contrasts(2) = coefTest(mdl, [0 -2 0 0 0 0 0 0 -2 0 0]);', valType{nv}));
%         eval(sprintf('stats.%s_quant_contrasts(3) = coefTest(mdl, [0 -2 0 0 0 0 0 0 0 -2 0]);', valType{nv}));
%         eval(sprintf('stats.%s_quant_contrasts(4) = coefTest(mdl, [0 -2 0 0 0 0 0 0 0 0 -2]);', valType{nv}));
% 
%         eval(sprintf('stats.%s_quant_contrasts(5) = coefTest(mdl, [0 -2 0 0 0 0 0 2 2 2 2]);', valType{nv}));



        

    end
end