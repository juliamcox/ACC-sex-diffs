function pairedEPSP 

basefilename = whereAreWe('figurecode');


try load(fullfile(basefilename,'figure_3','chr2_potential.mat'))
catch
    data = readtable(fullfile(basefilename, 'raw_data','figure_3', 'chr2_potentials.csv'));
    save(fullfile(basefilename,'figure_3','chr2_potentials.mat'),'data');
end

f = 'current ~ msn_type + (1|pair:mouse) + (1|mouse)';

data.msn_type = categorical(data.msn_type);
data.pair = categorical(data.pair);
data.mouse = categorical(data.mouse); 

lme = fitlme(data,f,'DummyVarCoding','effect');

