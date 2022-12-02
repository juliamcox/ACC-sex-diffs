function [cons, con_shift, time_back_orig, time_forward_orig, type1, shift_con, bsIDs,contVar,fullVar] = getEvents(regName,frameRate)

% cons: which events
% con_shift: will events be shifted to zero? (only for cases with 1 kernel size)
% time_back_orig: time before event in seconds (can be 1xnum events array or single value)
% time_forward_orig: time after event in seconds (can be 1x num events array or single value)
% type1: type of predictors 'spline' or 'time_shift' 
% shift_con: shift predictors? 0 or 1
% bsIDs: name of the basis set; format 'Xsec_Y' where X is the span in seconds and Y is the number of functions; can be an 1 X num events cell array

contVar = [];
fullVar = [];
% which events?
switch regName           

    case 'Basic3' %all events
        cons = {'NPTimes', 'LeverPresent','LeverTimesI', 'LeverTimesC', 'CSRew', 'CSNoRew'};
        con_shift = [0,1,0,0,1,1];
        %%%Set parameters        
        time_back_orig = 2; %seconds; time before event
        time_forward_orig = 6; %seconds; time after event
        type1 = 'spline'; %type of predictors convolved with basis set 'spline' or 'time_shift'
        shift_con = 1; %shift stimulus events so span starts at 0?
        bsIDs = cat(2,'8sec_25',['_' num2str(round(frameRate))]);       
    case 'outcome' %Figure 4E-F
        cons={'NPTimes','LeverPresent','LeverTimesI','LeverTimesC','CS','CSNoRew'};
        con_shift = [0,1,0,0,1,1];
        %%%Set parameters
        time_back_orig = 2; %seconds; time before event
        time_forward_orig = 6; %seconds; time after event
        type1 = 'spline'; %type of predictors convolved with basis set 'spline' or 'time_shift'
        shift_con = 1; %shift stimulus events so span starts at 0?
        bsIDs = cat(2,'8sec_25',['_' num2str(round(frameRate))]);
    case 'choice' %Figure 5C-D
        cons={'NPTimes','LeverPresent','LeverTimes','LeverTimesI','CSNoRew','CSRew'};
        con_shift = [0,1,0,0,1,1];
        %%%Set parameters
        time_back_orig = 2; %seconds; time before event
        time_forward_orig = 6; %seconds; time after event
        type1 = 'spline'; %type of predictors convolved with basis set 'spline' or 'time_shift'
        shift_con = 1; %shift stimulus events so span starts at 0?
        bsIDs = cat(2,'8sec_25',['_' num2str(round(frameRate))]);
    case 'stay' %Figure 4G-H
        cons={'NPTimes','LeverPresent','LeverTimesI','LeverTimesC','CSRewStay','CSRew','CSNoRewStay','CSNoRew'};
        con_shift = [0,1,0,0,1,1,1,1];
        %%%Set parameters
        time_back_orig = 2; %seconds; time before event
        time_forward_orig = 6; %seconds; time after event
        type1 = 'spline'; %type of predictors convolved with basis set 'spline' or 'time_shift'
        shift_con = 1; %shift stimulus events so span starts at 0?
        bsIDs = cat(2,'8sec_25',['_' num2str(round(frameRate))]);
        
 
end

