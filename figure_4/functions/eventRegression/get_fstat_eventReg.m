
function [F_vec,F_vec_noRefit,b,stats,mdl_reduced]=get_fstat_eventReg(X,y,preds_to_test_cell)

%%% Fig 2 %%%

% Caclulate f-statistic comparing full and reduced models with and without re-fitting and 5-fold cross validation 

% X: the matrix of predictors (without the constant term)
% y: the dependent variable
% preds_to_test_cell: a cell with the indexes of predictors to be tested, for example {[1 2],[3 5],6}


mdl = fitlm(X,y);

stats.full = mdl;
b = mdl.Coefficients.Estimate;


X = [ones(size(X,1),1) X]; % add intercept
n = size(X,1);             % number of observations
k = size(X,2);             % number of variables (including constant) for full model

%% Calculate F statistic with refitting of the reduced model 

sse_full = nansum((y - X * b).^2);
df_full = n - numel(b); 

for ne = 1:numel(preds_to_test_cell)
    % refit model without event's predictors 
    clear mdl 
    X_reduced = X;
    X_reduced(:,preds_to_test_cell{ne}+1) = [];
    
    mdl_reduced = fitlm(X_reduced(:,2:end),y);
    stats_reduced{ne} = mdl_reduced;
    b_reduced{ne} = mdl_reduced.Coefficients.Estimate;
    
    sse_reduced(ne) = nansum((y - X_reduced * b_reduced{ne}).^2);
    df_reduced = n - numel(b_reduced{ne}); 
    F_vec(ne) = ((sse_reduced(ne)-sse_full)/(df_reduced-df_full))/(sse_full/df_full);   
end

%stats.reduced = stats_reduced; 


%% Calculate F statistic without refitting the reduced model

X(isnan(y),:) = [];                  % remove NaNs
y(isnan(y))   = [];
u             = y - X * b;           % calculates residuals
s2            = u' * u / (n - k);    % estimate variance of error term (assuming homoskedasticity, independent observations) or sum(u.^2)/(n-k)
BCOV          = inv(X'*X) * s2;      % get covariance matrix of b assuming homoskedasticity of error term etc...
bse           = diag(BCOV).^.5;      % standard errors


for l=1:length(preds_to_test_cell)
    preds_to_test_cell{l} = preds_to_test_cell{l}+1; %because of the constant
    R = zeros(length(preds_to_test_cell{l}),k);
    for l2 = 1:length(preds_to_test_cell{l})
        R(l2, preds_to_test_cell{l}(l2))=1;
    end
    
    r = zeros(length(preds_to_test_cell{l}),1);          % Testing restriction: R * b = r
    num_restrictions = size(R, 1);
    F = (R*b - r)'*inv(R * BCOV * R')*(R*b - r) / num_restrictions;   % F-stat (see Hiyashi for reference)
    F_vec_noRefit(l) = F;
end




