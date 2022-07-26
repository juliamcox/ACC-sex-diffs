function tbl = statsTable(anovaIn, coeffIn, dfFlag, tflag,savename)


varNames = coeffIn(1:end,1);
if ~dfFlag
    if tflag
    tbl = table(varNames,'VariableNames',{'name'});
    est = cellfun(@(x,y) cat(2,num2str(x,'%0.3f'),char(177),num2str(y,'%0.2f')), coeffIn(2:end,2),coeffIn(2:end,3),'UniformOutput',false);
    est = cat(1,sprintf('\x03b2%sstandard error',char(177)),est);
    tbl.est = est;
    
    tstat = cellfun(@(x) num2str(x,'%0.2f'), coeffIn(2:end,4),'UniformOutput',false);
    tstat = cat(1,sprintf('F-stat(%d,%d)', 1,coeffIn{2,5}),tstat);
    tbl.tstat = tstat;
    
    pvals = [];
    for na = 2:numel(varNames)
        if coeffIn{na,6} < 0.001
            pvals = cat(1,pvals,{num2str(coeffIn{na,6},'%0.2e')});
        else
            pvals = cat(1,pvals,{num2str(coeffIn{na,6},'%0.2f')});
        end
    end
    pvals = cat(1,'p-value',pvals);
    tbl.pvals = pvals;
    
    ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,7),'UniformOutput',false);
    ci = cat(1,'lower',ci);
    tbl.ci1=ci;
    ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,8),'UniformOutput',false);
    ci = cat(1,'upper',ci);
    tbl.ci2=ci;
    else
    tbl = table(varNames,'VariableNames',{'name'});
    est = cellfun(@(x,y) cat(2,num2str(x,'%0.3f'),char(177),num2str(y,'%0.2f')), coeffIn(2:end,2),coeffIn(2:end,3),'UniformOutput',false);
    est = cat(1,sprintf('\x03b2%sstandard error',char(177)),est);
    tbl.est = est;
    
    fstat = cellfun(@(x) num2str(x,'%0.2f'), anovaIn(2:end,2),'UniformOutput',false);
    fstat = cat(1,sprintf('F-stat(%d,%d)', anovaIn{2,3},anovaIn{2,4}),fstat);
    tbl.fstat = fstat;
    
    pvals = [];
    for na = 2:numel(varNames)
        if anovaIn{na,5} < 0.001
            pvals = cat(1,pvals,{num2str(anovaIn{na,5},'%0.2e')});
        else
            pvals = cat(1,pvals,{num2str(anovaIn{na,5},'%0.2f')});
        end
    end
    pvals = cat(1,'p-value',pvals);
    tbl.pvals = pvals;
    
    ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,7),'UniformOutput',false);
    ci = cat(1,'lower',ci);
    tbl.ci1=ci;
    ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,8),'UniformOutput',false);
    ci = cat(1,'upper',ci);
    tbl.ci2=ci;
    end
else
    if tflag
        
        keyboard
    else
        
        tbl = table(varNames,'VariableNames',{'name'});
        est = cellfun(@(x,y) cat(2,num2str(x,'%0.3f'),char(177),num2str(y,'%0.2f')), coeffIn(2:end,2),coeffIn(2:end,3),'UniformOutput',false);
        est = cat(1,sprintf('\x03b2%sstandard error',char(177)),est);
        tbl.est = est;
        
        fstat = cellfun(@(x) num2str(x,'%0.2f'), anovaIn(2:end,2),'UniformOutput',false);
        fstat = cat(1,sprintf('F-stat(%d,%d)', anovaIn{2,3},anovaIn{2,4}),fstat);
        tbl.fstat = fstat;
        
        tbl.df1 = anovaIn(:,3);
        tbl.df2 = anovaIn(:,4);
        
        pvals = [];
        for na = 2:numel(varNames)
            if anovaIn{na,5} < 0.001
                pvals = cat(1,pvals,{num2str(anovaIn{na,5},'%0.3e')});
            else
                pvals = cat(1,pvals,{num2str(anovaIn{na,5},'%0.3f')});
            end
        end
        pvals = cat(1,'p-value',pvals);
        tbl.pvals = pvals;
        
        ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,7),'UniformOutput',false);
        ci = cat(1,'lower',ci);
        tbl.ci1=ci;
        ci = cellfun(@(x) num2str(x,'%0.2f'),coeffIn(2:end,8),'UniformOutput',false);
        ci = cat(1,'upper',ci);
        tbl.ci2=ci;
    end
    
end

writetable(tbl,savename)