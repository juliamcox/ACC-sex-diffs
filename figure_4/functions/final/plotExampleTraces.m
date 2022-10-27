function plotExampleTraces(recs)

fs=20;
rawFlag = 1;

basefilename = fullfile(whereAreWe('imaging'),recs);

load(fullfile(basefilename,'analyzed_neuron_final.mat'),'neuron')
try
    load(fullfile(basefilename,'analyzed_neuron_final.mat'),sprintf('dff_%d_raw%d',fs,rawFlag))
dff = eval(sprintf('dff_%d_raw%d',fs,rawFlag));
dff = nanzscore(dff,1);
catch
    dff = eval(sprintf('dff_%d',fs));
    load(fullfile(basefilename,'analyzed_neuron_final.mat'),sprintf('dff_%d',fs))
end
load(fullfile(basefilename,'Y.mat'),'rCrop','cCrop','Ysiz')

timeAxis = randi(size(dff,1)-(fs*45),1);
timeAxis = timeAxis:timeAxis+fs*45;
try
allNeurons = reshape(full(neuron.A),rCrop(2)-rCrop(1),cCrop(2)-cCrop(1),size(neuron.A,2));
catch
    allNeurons = reshape(full(neuron.A),size(neuron.Cn,1),size(neuron.Cn,2),size(neuron.A,2));
end
for nn = 1:size(dff,2)
    figure('Position', [680 441 307 537])
   subplot(2,1,1)
    plot(timeAxis,dff(timeAxis,nn), 'Color','k','LineWidth',1)
    subplot(2,1,2)
    try
    imagesc(reshape(neuron.A(:,nn),rCrop(2)-rCrop(1),cCrop(2)-cCrop(1)))
    catch
            imagesc(reshape(neuron.A(:,nn),size(neuron.Cn,1),size(neuron.Cn,2)));
    end
    figure(); imagesc(mean(allNeurons,3))
    pause()
    close all
end
