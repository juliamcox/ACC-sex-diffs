function [pval, chi2stat, posthoc] = chi2PropTest_multi(propTable,sigLevel)

% Perform Chi-square test to compare proportions in r x c prop table 

for nc = 1:size(propTable,2)
   for nr = 1:size(propTable,1)
      expected(nr,nc) = (sum(propTable(nr,:))*sum(propTable(:,nc)))/sum(sum(propTable));  
   end
end

chi2stat = sum(sum(((propTable-expected).^2)./expected));
df = (size(propTable,1)-1)*(size(propTable,2)-1);  

pval = 1-chi2cdf(chi2stat,df); 

% Perform pairwise comparisons with the Marascuillo procedure

props = reshape(propTable./sum(propTable,2),1,size(propTable,1)*size(propTable,2));
nsample = reshape(repmat(sum(propTable,2),1,size(propTable,2)),1,size(propTable,1)*size(propTable,2));

counter = 1; 
for n = 1:numel(props)
   for nn = 1+n:numel(props) 
       critRange(counter) = sqrt(chi2inv(1-sigLevel,df))*sqrt(((props(n)*(1-props(n)))/nsample(n)) + ((props(nn)*(1-props(nn)))/nsample(nn)));
       propDiff(counter) = abs(props(n)-props(nn)); 
       comparison(counter,:) = [n nn]; 
       counter = counter+1;
   end
end

sigtable = abs(propDiff)>critRange;

posthoc.comparison = comparison;
posthoc.propDiff = propDiff;
posthoc.critValue = critRange;
posthoc.df = df;
posthoc.sigtable = sigtable;
posthoc.props = props;