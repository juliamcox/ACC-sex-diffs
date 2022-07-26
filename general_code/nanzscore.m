function data = nanzscore(data,dim)

if nargin < 2 
    [~,dim] = (max(size(data)));
end


mu = nanmean(data,dim);
sd = nanstd(data,[],dim);

data = (data - mu)./sd; 