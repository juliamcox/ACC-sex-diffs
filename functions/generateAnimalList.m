function [aids] = generateAnimalList(cohort)

% Animal IDs for various cohorts

switch cohort
    case 'estrous'
        aids = {'A10';'A11';'A12';'A13';'A14';'A20';'A21';'A22';'A23';'A24';'A31';'A32';'A33'};
        % not cycling 'A30';
    case 'ACC_DMS_nphr_male'
        aids = {'N1';'N2';'N3';'T46';'106';'107';'108';'109';'110'};
    case 'ACC_DMS_nphr_female'
        aids = {'N6';'N7';'N8';'T52';'T96';'T97';'T98';'T99'}; 
    case 'ACC_DMS_nphr_yfp_male'
        aids= {'N4';'N5';'T49';'T48';'104';'105'};
    case 'ACC_DMS_nphr_yfp_female'
        aids ={'N9';'N10';'T41';'T44';'T50';'T54'};
    
    case 'DMS_nphr_d1_male'
        aids = {'D16';'D15';'D23';'D32';'D28'};
    case 'DMS_nphr_d1_female'
        aids = {'D9';'D10';'D18';'D37';'D38';'D39'};
    case 'DMS_nphr_d2_male'
        aids = {'I15';'I14';'I19';'I17';'I13';'I24';'I30'};
    case 'DMS_nphr_d2_female'
        aids = {'I39';'I40';'I41';'I47';'I48';'I49';'I42';'I44'};
        
    case 'DMS_yfp_female'
        aids = {'I45';'I46';'I50';'I53';'I51';'I52';'D17'};
    case 'DMS_yfp_male'
        aids = {'D20';'D29';'D31';'D33';'I18';'I27';'I28';'I29';'I33';'I31'};
        
    
    case 'imaging_male' % names for behavior analyses
        aids = {'T1';'T3';'T12';...
            'T14';'T60';'T63';'102';'103';'101'};
    case 'ACC_DMS_imaging_male'
        %         % 'T72';'recording_20191202_112218'; performance too bad for value model
        aids = {'T1';'T3';'T12';...
            'T14';'T60';'T63';'102';'103';'101'};
        recs = {'recording_20161204_142407';'recording_20161130_173418';'recording_20170310_125113';...
            'recording_20170220_115346';'recording_20190725_145401';'recording_20190627_122511';'recording_20200721_140825';'recording_20200711_135354';'recording_20200722_152930'};
        temp  = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
        % T60 'recording_20190625_132507'recording_20190725_145401
        
    case 'ACC_DMS_imaging_male_2'
        aids = {'102';'103';'101';
            'T1';'T3';'T12';...
            'T14';'T60';'T63'};
        recs = {'recording_20200709_145812';'recording_20200723_121031';'recording_20200718_143237'};
        temp  = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
        
    case 'imaging_female'
        aids = {'T6';'T19';'T57';'T55';'T95';'T94'};
        
    case 'ACC_DMS_imaging_female'
 
        aids = {'T6';'T19';'T57';'T55';'T95';'T94'};
        recs = {'recording_20161207_154909' ;'recording_20170910_135133';'recording_20190731_123853';'recording_20190627_162644';'recording_20200722_141305';'recording_20200717_160625'}; 
        %Alternate T6 recording is  'recording_20161211_120749'temp  = [];
        %alternate T57 ' although there seems
        %torecording_20190730_144234' 
        %be a problem with the sync file
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_female_2'
        
        aids = {
            'T6';'T19';'T57';'T55';'T95';'T94'};
        recs = {};
        %Alternate T6 recording is  'recording_20161211_120749'temp  = [];
        %alternate T57 ' although there seems
        %torecording_20190730_144234'
        %be a problem with the sync file
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    
    case 'imaging_female_group'
        aids = {'T65';'T68';'T67'};
    case 'imaging_male_group'
        aids = {'T71';'T73';'T74'};
    case 'ACC_DMS_imaging_female_group'
        aids = {'T65';'T68';'T67'};
        recs = {'recording_20191118_123924';'recording_20191106_123959';'recording_20191202_123714'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_male_group'
        aids = {'T71';'T73';'T74'};
        recs = {'recording_20191106_112952';'recording_20191120_101441';'recording_20191120_112814'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_female_shock_group'
        aids = {'T68';'T65';'T67'};
        recs = {'recording_20191107_130148';'recording_20191211_155555';'recording_20191213_124253'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_male_shock_group'
        aids = {'T71';'T73'};
        recs = {'recording_20191107_135502';'recording_20200106_113450'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_female_shock_group_day2'
        aids = {'T68';'T65';'T67'};
        recs = {'recording_20191108_112816';'recording_20191213_144934';'recording_20191214_150012';};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_male_shock_group_day2'
        aids = {'T71';'T73'};
        recs = {'recording_20191108_121936';'recording_20200106_113450'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
    case 'ACC_DMS_imaging_male_shock'
        aids = {'T72'};
        recs = {'recording_20191213_180641'};
        temp = [];
        for na = 1:numel(aids)
            temp = cat(1,temp, {fullfile(aids{na},recs{na})});
        end
        aids = temp;
end



