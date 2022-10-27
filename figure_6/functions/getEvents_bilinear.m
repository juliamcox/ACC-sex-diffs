function [events, kernelType, gains, eventDur,nbasis,numConsec] = getEvents_bilinear(regName)

%events: cell array of names of events:
%'nosePokeEntry';'leverPresentation';rewardEntry;CSMinus;'outcome';'leverPress';'leverPressIpsi';'leverPressContra'
%'leverPress
%gains: cell array of trial parameters to use to fit gains on kernels
%kernelType: which type of basis set, if any, to use  'bspline','cosine' or 'timeLag'
%eventDur: time before and after event to include in kernel, in seconds
%nbasis: number of basis sets per event if kernelType is bspline or cosine, otherwise empty
%numConsec: number of consecutive rewards or no rewards (if applicable, and only for plotting)
numConsec = 0; 
% which events?
switch regName           
   
    case 'v10'
        events     = {'nosePokeEntry','leverPresentation','leverPressIpsi','leverPressContra','outcome'};
        kernelType = 'bspline';
        gains      = {{'qDiff2'}; {'qDiff2'};{'qDiff2'};{'qDiff2'};{'qDiff2';'qDiff2_next'}};
        eventDur   = [2 6; 0 8; 2 6; 2 6; 0 8];
        nbasis     = [25;25;25;25;25];        
    case 'v26'
        events     = {'nosePokeEntry','leverPresentation','leverPress','CSPlus','CSMinus'};
        kernelType = 'bspline';
        gains      = {{'qDiff'}; {'qDiff'};{'qDiff'};{'qDiff';'qDiff_next'};{'qDiff';'qDiff_next'}};
        eventDur   = [2 6; 0 8; 2 6; 0 8; 0 8];
        nbasis     = [25;25;25;25;25;25];
    case 'v12'
        events     = {'nosePokeEntry','leverPresentation','leverPressIpsi','leverPressContra','outcome'};
        kernelType = 'bspline';
        gains      = {{'qTot'}; {'qTot'};{'qTot'};{'qTot'};{'qTot'; 'qTot_next'}};
        eventDur   = [2 6; 0 8; 2 6; 2 6; 0 8];
        nbasis     = [25;25;25;25;25];
end