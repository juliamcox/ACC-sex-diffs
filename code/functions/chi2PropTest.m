function [pval, chi2stat] = chi2PropTest(k1,k2,n1,n2,correctFlag)

% Perform Chi-square test to compare 2 proportions
% k1 number of occurrences for group 1
% k2 number of occurrences for group 2
% n1 total size of group 1
% n2 total size of group 2 

p1 = k1/n1; %proportion of success in group 1
p2 = k2/n2; %proportion of success in group 2


% check that there are enough samples in each comparison group
if min([n1*p1,n1*(1-p1),n2*p2,n2*(1-p2)]) < 5
    disp('Error: fewer than 5 samples in at least one group')
end

df = 1; 

p = (k1+k2)/(n1+n2); % total proportion 

% expected counts
k01 = n1*p;
k02 = n2*p;

observed = [k1 n1-k1 k2 n2-k2]; %success group 1, failure group 1, success group 2, failure group 2
expected = [k01 n1-k01 k02 n2-k02];
if correctFlag
    chi2stat = sum((abs(observed-expected)-.5).^2 ./expected); %Yates continuity correction
else
    chi2stat = sum(((observed-expected).^2)./expected);
end
pval = 1-chi2cdf(chi2stat,df); 