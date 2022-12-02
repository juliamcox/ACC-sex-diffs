function stats=ephysVariability(data)


pairs_m = unique(data.pair(data.sex==1));
pairs_f = unique(data.pair(data.sex==2));


for np = 1:numel(pairs_m)
    diff_m(np) = data.current(data.pair==pairs_m(np)&data.msn_type==2)-data.current(data.pair==pairs_m(np)&data.msn_type==1);
end

for np = 1:numel(pairs_f)
    diff_f(np) = data.current(data.pair==pairs_f(np)&data.msn_type==2)-data.current(data.pair==pairs_f(np)&data.msn_type==1);
end

stats.var_m = var(diff_m);
stats.var_f = var(diff_f);

[stats.h,stats.p,~,stats.stats] = vartest2(diff_f,diff_m);