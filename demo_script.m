%% Demo script for phase amplitude coupling toolbox.
%
%
%
% This toolbox supports the following methods :
% 
% * Mean Vector Length Modulation Index _Canolty et al. 2006_
% * Kullback-Leibler Modulation Index _Tort et al. 2010_
% * Generalized Linear Model _Penny et al. 2008_
% * MIPAC _Martinez-Cancino et al. 2018_
% 
%
%
% Authors: Ramon Martinez Cancino,SCCN/INC, UCSD 2018
%
% Copyright (C) 2018 Ramon Martinez Cancino,UCSD, INC, SCCN

%% Load simulated single trial PAC data
% A single trial PAC simulated signal with phase at 8 Hz coupled with the amplitude ate 60Hz 
% (see Martinez-Cancino et al., 2018 for details on the simulation).
% Duration of the signals is 7s at a sampling rate of 500Hz. The signals
% was splitt

EEG = pop_loadset('filename','strial_sim_pac_ph8amp60.set','filepath',pwd);

% Visualizing the simulated data
eegplot(EEG.data,'srate', 500,'winlength',7)

%% Compute PAC using Instanteaneous MIPAC method
% Parameters
mipacvarthreshval = 0.05; % Variance Threshold for MIPAC iterations
phaserange        = [5 12];
nfreqsphase       = 8;
amprange          = [30 100];
nfreqsamp         = 15;

[pacval, ~, ~, ~, ~, ~,~, instmipac_pacstruct] = eeg_pac(EEG.data(1,:)', EEG.data(1,:)', EEG.srate,...
                                                         'methodpac',      'instmipac',...
                                                         'freqs',          phaserange,...                                                        
                                                         'nfreqs1',        nfreqsphase,...
                                                         'freqs2',         amprange,... 
                                                         'nfreqs2',        nfreqsamp,...                                                                                                             
                                                         'timefreq',       0, ...
                                                         'butterorder',    6, ...
                                                         'mipacvarthresh', mipacvarthreshval,...
                                                         'alpha',          []);
                                                     
% Visualizing MIPAC Time series (tool for visualization under development)
figure('units', 'normalized','position',[0.1203 0.3944 0.5797 0.2722]); 
instmipac_pacstruct.params.freqs_amp
instmipac_pacstruct.params.freqs_phase
plot(instmipac_pacstruct.instmipac.times,squeeze(instmipac_pacstruct.instmipac.pacval(5,8,:)));
xlim(minmax(instmipac_pacstruct.instmipac.times)); xlabel('Time (s)');

% Visualizing MIMI comodulogram
figure('units', 'normalized','position',[0.1187 0.3972 0.2484 0.4083]); 
imagesc(mvlmi_pacstruct.params.freqs_phase,mvlmi_pacstruct.params.freqs_amp,flipud(mean(instmipac_pacstruct.instmipac.pacval,3)'));
set(gca,'YDir','normal');
%% MVLMI (no significance analysis)
% Parameters
phaserange        = [5 12];
nfreqsphase       = 8;
amprange          = [30 100];
nfreqsamp         = 15;

[pacval, ~, ~, ~, ~, ~,~, mvlmi_pacstruct] = eeg_pac(EEG.data(1,:)', EEG.data(1,:)', EEG.srate,...
                                                         'methodpac', 'mvlmi',...
                                                         'freqs',     phaserange ,...
                                                         'nfreqs1',   nfreqsphase,...
                                                         'freqs2',    amprange,...
                                                         'nfreqs2',   nfreqsamp,...
                                                         'timefreq',  0,...
                                                         'alpha',     []);
                                                     
% Visualizing (tool for visualization under development)
figure('units', 'normalized','position',[0.1187 0.3972 0.2484 0.4083]); 
imagesc( mvlmi_pacstruct.params.freqs_phase, mvlmi_pacstruct.params.freqs_amp,flipud(mvlmi_pacstruct.mvlmi.pacval'));
set(gca,'YDir','normal');                                                     
%% GLM (no significance analysis)
% Parameters
phaserange        = [5 12];
nfreqsphase       = 8;
amprange          = [30 100];
nfreqsamp         = 15;

[pacval, ~, ~, ~, ~, ~,~, glm_pacstruct] = eeg_pac(EEG.data(1,:)', EEG.data(1,:)', EEG.srate,...
                                                         'methodpac', 'glm',...
                                                         'freqs',     phaserange ,...
                                                         'nfreqs1',   nfreqsphase,...
                                                         'freqs2',    amprange,...
                                                         'nfreqs2',   nfreqsamp,...
                                                         'timefreq',  0,...
                                                         'alpha',     []);    
                                                     
% Visualizing (tool for visualization under development)
figure('units', 'normalized','position',[0.1187 0.3972 0.2484 0.4083]);
imagesc(glm_pacstruct.params.freqs_phase, glm_pacstruct.params.freqs_amp,flipud(glm_pacstruct.glm.pacval'));   
set(gca,'YDir','normal');
%% KLMI (no significance analysis)
% Parameters
phaserange        = [5 12];
nfreqsphase       = 8;
amprange          = [30 100];
nfreqsamp         = 15;

[pacval, ~, ~, ~, ~, ~,~, klmi_pacstruct] = eeg_pac(EEG.data(1,:)', EEG.data(1,:)', EEG.srate,...
                                                         'methodpac', 'klmi',...
                                                         'freqs',     phaserange ,...
                                                         'nfreqs1',   nfreqsphase,...
                                                         'freqs2',    amprange,...
                                                         'nfreqs2',   nfreqsamp,...
                                                         'timefreq',  0,...
                                                         'alpha',     []);   
                                                     
% Visualizing (tool for visualization under development)
figure('units', 'normalized','position',[0.1187 0.3972 0.2484 0.4083]); 
imagesc(klmi_pacstruct.params.freqs_phase, klmi_pacstruct.params.freqs_amp,flipud(klmi_pacstruct.klmi.pacval'));
set(gca,'YDir','normal');
%==========================================================================
%==========================================================================