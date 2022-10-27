function  activeIdx = findActiveNeurons(recs,fbasename)

 

load(fullfile(fbasename,recs,'analyzed_neuron_final.mat'), 'neuron'); 
info = matfile(fullfile(fbasename,recs,'info.mat'));
frameRate = info.frameRate; 

spikeRate = mean(neuron.S,2) .* frameRate .* 60; 

activeIdx = find(spikeRate>1); 

save(fullfile(fbasename,recs,'activeNeurons.mat'),'activeIdx'); 
