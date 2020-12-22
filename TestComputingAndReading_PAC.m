% std_pacprecomp  and std_readpacdata testing cases.

% Loading STUDY.
% ERSP for channels and components are already computed. Clusters are
% computed as well.
% Here testing with a two-subject STUDY from the 5-subject EEGLAB testing STUDY 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[STUDY ALLEEG] = pop_loadstudy('filename', 'n400clustedit.study', 'filepath', '/Users/amon-ra/WORK/Data/toplay/studytest');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
eeglab redraw


%% Generating and reading PAC measures computed with ERMIPAC
%---------------------------------------------------------
%---------------------------------------------------------
%% Calling std_precomp channels [DO NOT RUN. TAKES A CONSIDERABLE AMOUNT OF TIME. FIX NEEDED]
chanorcomp = 1;
chanindx = {'Fp1' , 'F3' , 'F5'};
freqphase = [3 6 9];
freqamp   = [31 45 70];
pacparams = {'method', 'ermipac'};
STUDYout = std_pacprecomp(STUDY, ALLEEG,freqphase,freqamp, chanorcomp,chanindx, 'method', 'ermipac', 'pacparams',{'k', 5});

%% Calling std_precomp IC clusters [DO NOT RUN. TAKES A CONSIDERABLE AMOUNT OF TIME. FIX NEEDED]
chanorcomp = 2;
chanindx = {'Cls 10'};
freqphase = [3 6 9];
freqamp   = [31 45 70];
pacparams = {'method', 'ermipac'};
STUDYout = std_pacprecomp(STUDY, ALLEEG,freqphase,freqamp, chanorcomp,chanindx, 'method', 'ermipac', 'pacparams',{'k', 5});


%% Test std_readpacdata
%------------------------
%% Reading channel data (single channel all subjects)
[STUDY, pac1, alltimes1, freq1_1, freq2_1, ~, fileparams] = std_readpacdata(STUDY, ALLEEG, 'channels', {'Fp1'}, 'timerange',[],'subject', '', 'singletrials', 'on', 'design', 1, 'datatype', 'MIPAC');
% Output here must be:  2×1 cell array
%                                       {3×3×17×458 single}
%                                       {3×3×17×463 single}

%% Reading channel data (single channel single subject)
[STUDY, pac1, alltimes1, freq1_1, freq2_1, ~, fileparams] = std_readpacdata(STUDY, ALLEEG, 'channels', {'Fp1'}, 'timerange',[],'subject', 'S02', 'singletrials', 'on', 'design', 1, 'datatype', 'MIPAC');
% Output here must be:  2×1 cell array
%                                       {3×3×17×223 single}
%                                       {3×3×17×231 single}


%% Reading channel data (Multiple channel) 
[STUDY, pac2, alltimes2, freq1_2, freq2_2, ~, fileparams] = std_readpacdata(STUDY, ALLEEG, 'channels', {'Fp1' , 'F3' , 'F5', 'F7'}, 'timerange',[],'subject', '', 'singletrials', 'on', 'design', 1, 'datatype', 'MIPAC');
% Output here must be:  2×1 cell array
%                                       {3×3×17x4x458 single}
%                                       {3×3×17×4x463 single}

%% Reading component data (single cluster all ICs) 
[STUDY, pac3, alltimes3, freq1_3, freq2_3, ~, fileparams] = std_readpacdata(STUDY, ALLEEG, 'clusters', 10, 'timerange',[],'component', '', 'singletrials', 'on', 'design', 1, 'datatype', 'MIPAC');
% Output here must be:  2×1 cell array
%                                     {3×3×17×1585 single}   
%                                     {3×3×17×1619 single}   

%% Reading component data (single cluster single IC) 
[STUDY, pac3, alltimes3, freq1_3, freq2_3, ~, fileparams] = std_readpacdata(STUDY, ALLEEG, 'clusters', 10, 'timerange',[],'component',1, 'singletrials', 'on', 'design', 1, 'datatype', 'MIPAC');
% Output here must be:  2×1 cell array
%                                     {3×3×17×1585 single}   
%                                     {3×3×17×1619 single}   

%% Generating and reading PAC measures computed with KLMI (representative of common PAC measures)
%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------

%% Calling std_precomp channels
chanorcomp = 1;
chanindx = {'Fp1' , 'F3' , 'F5', 'F7'};
freqphase = [3 6 9];
freqamp   = [31 45 70];
pacparams = {'method', 'klmi'};
STUDY = std_pacprecomp(STUDY, ALLEEG,freqphase,freqamp, chanorcomp,chanindx, 'method', 'klmi');

%% Calling std_precomp IC clusters
chanorcomp = 2;
chanindx = {'Cls 10', 'Cls 11'};
freqphase = [3 6 9];
freqamp   = [31 45 70];
pacparams = {'method', 'klmi'};
STUDY = std_pacprecomp(STUDY, ALLEEG,freqphase,freqamp, chanorcomp,chanindx, 'method', 'klmi');

%% Test std_readpacdata
%------------------------
%% Reading channel data (single channel)
[STUDY, pac1, alltimes1, freq1_1, freq2_1, ~, fileparams] = std_readpacdata(STUDYout, ALLEEG, 'channels', {'Fp1'}, 'timerange',[],'subject', '', 'singletrials', 'off', 'design', 1, 'datatype', 'PAC');
% Output here must be:  2×1 cell array
%                                       {3×3×20x2 single}
%                                       {3×3×20x2 single}
%% Reading channel data (Multiple channel) 
[STUDY, pac2, alltimes2, freq1_2, freq2_2, ~, fileparams] = std_readpacdata(STUDYout, ALLEEG, 'channels', {'Fp1' , 'F3' , 'F5', 'F7'}, 'timerange',[],'subject', '', 'singletrials', 'off', 'design', 1, 'datatype', 'PAC');
% Output here must be:  2×1 cell array
%                                       {3×3×20×4x2 single} 
%                                       {3×3×20×4x2 single}  

%% Reading component data (single cluster) 
[STUDY, pac3, alltimes3, freq1_3, freq2_3, ~, fileparams] = std_readpacdata(STUDYout, ALLEEG, 'clusters', 10, 'timerange',[],'subject', '', 'singletrials', 'off', 'design', 1, 'datatype', 'PAC');
% Output here must be:  2×1 cell array
%                                     {3×3×20×7 single} 
%                                     {3×3×20×7 single} 
 