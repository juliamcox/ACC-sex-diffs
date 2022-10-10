idx = strfind(neuron.P.log_file,'/');
if isempty(idx)
    idx = strfind(neuron.P.log_file, '\');
end
neuron.P.log_file=fullfile(cd,neuron.P.log_file(idx(end-3):end));
neuron.P.mat_file= matfile((fullfile(cd,'data_source_extraction/data_64_64_18.mat')));
idx = strfind(neuron.P.log_data,'/');
if isempty(idx)
    idx = strfind(neuron.P.log_data, '\');
end
neuron.P.log_data= fullfile(cd,neuron.P.log_data(idx(end-3):end));