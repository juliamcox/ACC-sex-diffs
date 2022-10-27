function data = concatValue(data,qLearn,aids)

basefilename = fullfile(whereAreWe('bucket'), 'Operant');

% Create array of Q-values for the sessions in data.daylist
try
    flist = {qLearn.fList.name};
catch
    flist = qLearn.fList;
end
for nd = 1:numel(data.daylist)
    startIdx = strfind(data.daylist{nd}, '\');
    if isempty(startIdx)
        startIdx = strfind(data.daylist{nd},'/');
    end
    endIdx = strfind(data.daylist{nd},'_');
    if isempty(startIdx)
        fIdx = find(contains(flist,data.daylist{nd}(1:endIdx(end)-1)));
    else
        fIdx = find(contains(flist,data.daylist{nd}(startIdx(end)+1:endIdx(end)-1)));
    end
    if numel(fIdx) > 1
        keyboard
    end
    if isempty(fIdx)
        % If the session didn't get run with the Qlearning model 
        if ~isempty(startIdx)
            warning(sprintf('No Q-values for session %s', data.daylist{nd}(startIdx(end)+1:endIdx(end)-1)));
            data2 = load(fullfile(basefilename,aids, data.daylist{nd}(startIdx(end)+1:end)));
        else
            warning(sprintf('No Q-values for session %s', data.daylist{nd}(1:endIdx(end)-1)));
            data2 = load(fullfile(basefilename,aids, data.daylist{nd}(1:end)));
        end
        
        QRight_temp{nd} = nan(size(data2.data.choice'));
        QLeft_temp{nd} = nan(size(data2.data.choice'));
        QIpsi_temp{nd} = nan(size(data2.data.choice'));
        QContra_temp{nd} = nan(size(data2.data.choice'));
        probIpsi_temp{nd} = nan(size(data2.data.choice'));
        probContra_temp{nd} = nan(size(data2.data.choice'));
    else
        trialIdx = find(qLearn.session == fIdx);
        if isempty(startIdx)
            data2 = load(fullfile(basefilename,aids, data.daylist{nd}(1:end)));
        else
            data2 = load(fullfile(basefilename,aids, data.daylist{nd}(startIdx(end)+1:end)));
        end
        try
        QRight_temp{nd} = qLearn.QRight(trialIdx);
        QLeft_temp{nd} = qLearn.QLeft(trialIdx);
        
        if contains(data2.data.laserSide,'L')
            QIpsi_temp{nd} = QLeft_temp{nd};
            QContra_temp{nd} = QRight_temp{nd};
            probIpsi_temp{nd} = qLearn.probLeft(trialIdx);
            probContra_temp{nd} = qLearn.probRight(trialIdx); 
        elseif contains(data2.data.laserSide,'R')
            QIpsi_temp{nd} = QRight_temp{nd};
            QContra_temp{nd} = QLeft_temp{nd};
            probIpsi_temp{nd} = qLearn.probRight(trialIdx);
            probContra_temp{nd} = qLearn.probLeft(trialIdx); 
        else
            keyboard
        end
        
        catch
            keyboard
        end
     
    end
end

% Create vector of Q-values
QRight = nan(size(data.choice));
QLeft = nan(size(data.choice));
QIpsi = nan(size(data.choice));
QContra = nan(size(data.choice)); 
QRight_temp = cell2mat(QRight_temp');
QLeft_temp = cell2mat(QLeft_temp');
QIpsi_temp = cell2mat(QIpsi_temp');
QContra_temp = cell2mat(QContra_temp'); 
probContra = nan(size(data.choice));
probIpsi = nan(size(data.choice));
probIpsi_temp = cell2mat(probIpsi_temp');
probContra_temp = cell2mat(probContra_temp'); 

try
thisIdx = find(data.choice~=-10); 
QRight(thisIdx) = QRight_temp;
catch
    keyboard
end
QLeft(thisIdx) = QLeft_temp;
QContra(thisIdx) = QContra_temp;
QIpsi(thisIdx) = QIpsi_temp; 

data.probIpsi(thisIdx) = probIpsi_temp;
data.probContra(thisIdx) = probContra_temp; 

data.QRight = QRight;
data.QLeft = QLeft; 
data.QContra = QContra;
data.QIpsi = QIpsi; 