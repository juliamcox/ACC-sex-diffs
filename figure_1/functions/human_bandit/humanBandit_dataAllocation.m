function humanBandit_dataAllocation(cohort)

basefilename = fullfile(whereAreWe('figurecode'),'raw_data','figure_1','human_bandit');
basefilename2 = fullfile(whereAreWe('figurecode'),'processed_data','human_bandit');
if ~isfolder(fullfile(basefilename, sprintf('%s_f',cohort)))
    mkdir(fullfile(basefilename, sprintf('%s_f',cohort)))
    mkdir(fullfile(basefilename, sprintf('%s_m',cohort)))
end

%find file names 
cd(fullfile(basefilename,cohort));
fnames = dir('*.json');
fnames = {fnames(:).name};

%convert to mat files 
dataConvert_humanBandit(cohort);

%move files to appropriate folder
status=[];
stats2=[];
noReportCounter = 0;
otherCounter = 0;
for nf = 1:numel(fnames)
        thisFname = sprintf('%s.mat',fnames{nf}(1:end-5));
        load(fullfile(basefilename,cohort,thisFname));
        temp = data.raw{end};
        if contains(lower(temp.responses.gender_categorical),'female')
            status=cat(1,status,movefile(fullfile(basefilename,cohort,fnames{nf}),fullfile(basefilename,sprintf('%s_f',cohort))));
            status2=cat(1,status,movefile(fullfile(basefilename2,cohort,thisFname),fullfile(basefilename2,sprintf('%s_f',cohort))));
        elseif contains(lower(temp.responses.gender_categorical),'male')
            status=cat(1,status,movefile(fullfile(basefilename,cohort,fnames{nf}),fullfile(basefilename,sprintf('%s_m',cohort))));
            status2=cat(1,status,movefile(fullfile(basefilename2,cohort,thisFname),fullfile(basefilename2,sprintf('%s_m',cohort))));
        elseif contains(lower(temp.responses.gender_categorical),'Rather not say')
            noReportCounter = noReportCounter+1;
        elseif contains(lower(temp.responses.gender_categorical),'Other')
            otherCounter = otherCounter+1;
        else
            noReportCounter = noReportCounter+1;
            keyboard
        end
end

if sum(status)~=numel(status) || sum(status2)~=numel(status2)
    keyboard
end

fprintf('no report: %d',noReportCounter)
fprintf('other report: %d',otherCounter)
    