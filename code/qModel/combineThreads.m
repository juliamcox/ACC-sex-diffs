fbasename = fullfile(whereAreWe('figurecode'),'processed_data','qModel_all');

flist = {'run1_2';'run3';'run4'};

tempSess = [];
tempAnimal = [];
for nf = 1:numel(flist)
    load(fullfile(fbasename,flist{nf},'alphas.mat'));
    load(fullfile(fbasename,flist{nf},'alphas_mice.mat'));
    tempSess = cat(1,tempSess,vect);
    tempAnimal = cat(1,tempAnimal,vect);
end
alphas = tempSess;
save(fullfile(fbasename,'alphas.mat'),'alphas');
alphas = tempAnimal;
save(fullfile(fbasename,'alphas_mice.mat'),'alphas');

tempSess = [];
tempAnimal = [];
for nf = 1:numel(flist)
    load(fullfile(fbasename,flist{nf},'sides.mat'));
    load(fullfile(fbasename,flist{nf},'sides_mice.mat'));
    tempSess = cat(1,tempSess,vect);
    tempAnimal = cat(1,tempAnimal,vect);
end
sides = tempSess;
save(fullfile(fbasename,'sides.mat'),'sides');
sides = tempAnimal;
save(fullfile(fbasename,'sides_mice.mat'),'sides');

tempSess = [];
tempAnimal = [];
for nf = 1:numel(flist)
    load(fullfile(fbasename,flist{nf},'stays.mat'));
    load(fullfile(fbasename,flist{nf},'stays_mice.mat'));
    tempSess = cat(1,tempSess,vect);
    tempAnimal = cat(1,tempAnimal,vect);
end
stays = tempSess;
save(fullfile(fbasename,'stays.mat'),'stays');
stays = tempAnimal;
save(fullfile(fbasename,'stays_mice.mat'),'stays');

tempSess = [];
tempAnimal = [];
for nf = 1:numel(flist)
    load(fullfile(fbasename,flist{nf},'betas.mat'));
    load(fullfile(fbasename,flist{nf},'betas_mice.mat'));
    tempSess = cat(1,tempSess,vect);
    tempAnimal = cat(1,tempAnimal,vect);
end
betas = tempSess;
save(fullfile(fbasename,'betas.mat'),'betas');
betas = tempAnimal;
save(fullfile(fbasename,'betas_mice.mat'),'betas');