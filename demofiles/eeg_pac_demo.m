% Demo script for phase amplitude coupling toolbox.
%
% Generates artificial data then applies methods to detect phase amplitude 
% coupling.
% Can generate multiple trial on single trial data. In order for the
% methods to work on single trial data it is necessary to add noise to the
% signal (snrval ~= Inf). For multiple trials we recommend adding noise and 
% shifting slightly the data (maxshift > 0). This script will give plots
% of times of modulation and no modulation.
% Note : alpha = [] will not compute the statistics of the pac and decrease
% significantly the computation time.
%
%
% Author: Joseph Heng and Ramon Martinez Cancino, EPFL, SCCN/INC, UCSD 2016
%
% Copyright (C) 2002  Joseph Heng and Ramon Martinez Cancino, EPFL, UCSD, 
% INC, SCCN
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% Parameters
% These parameters can be modified.

% Signal parameters
fc            = 60;   % Frequency of the carrier wave. Recommended: [60]
fm            = 8;    % Frequency of the modulating wave. Recommended: [60]
max_time      = 10;    % Time of the simulation. Recommended: single trial [6], multiple trials [2]
n_trials      = 200;    % Number of trials. Recommended: single trial [1], multiple trials [1000]
s_rate        = 500;  % Sampling rate. Recommended: single trial [500], multiple trials [250]
padtime       = 0;  % Padding time at the start and end of the signal. 
                      % Should increase with sampling size.
                      % Recommended: 0.5
snrval        = 5;    % Signal to noise ratio. Recommended: [5]
maxshift      = 10;   % Maximum shifts (jitter) injected into the data.Recommended: single trial [1], multiple trials [10]
nsegm         = 5;    % Number of segments in the data. Each segment is a block with our without modulation


% Method parameters                      
nfreqsphase = 10;       % Number of frequencies used for the phase. Recommended: [10]
nfreqsamp   = 8;        % Number of frequencies used for the amplitude. Recommended: [10]
alpha       = [];       % Alpha value used for statistics. If alpha = [] the statistics are not computed.
                        % Recommended: [0.05] (or [] for less computationtime)
compute_mvl = true;     % Compute the Mean Vector Length Modulation Index
compute_kl  = true;     % Compute the Kullback Leibler Modulation Index
compute_glm = true;     % Compute the GLM Modulation Index
phaserange  = [4 25];   % Range of phases to be used
                        % Recommended: [4 25]
amprange    = [30 100]; % Range of amplitudes to be used  Recommended: [30 100]
nsurrogate  = 200;       % Number of surrogates. Recommended: [200].
nbinskl     = 18;       % Number of bins for the Kullback Leibler Modulation
                        % index. Recommended: [8]

% Do not modify the code beyond this point
%% Generate data
clear data_pac
tlimits      = [0 max_time];

s = randi(maxshift, n_trials,1);

for i=1:n_trials
     plot_flag = i == 1;
     [datatmp,t_outtmp]  = generate_pac_signal(fc,fm,tlimits,'Ac',5,'Am',1,'fs',s_rate,'cpfunc','block','blockamp',1, 'nsegm',nsegm,'plot_flag', plot_flag,'padtime',padtime,'snr',snrval,'m',1);
     datatmp = [datatmp(s(i):end) datatmp(1:s(i)-1)];
     data_pac(i,:) = [datatmp(end-ceil(maxshift/2):end) datatmp(1:end-ceil(maxshift/2)-1)]; 
end

%data_pac = data_pac(length(data_pac)/2:end);
%data_pac = data_pac(1:length(data_pac)/2);

timepointpac   = padtime + 3 * max_time/nsegm/2;
timepointnopac = padtime +     max_time/nsegm/2;
  
%% Compute PAC MVLMI
if compute_mvl,
    method = 'mvlmi';
    tic;
    [mvlmipacval, ~, ~, ~, ~, ~,~,~, mvlmipacstruct] = eeg_pac(data_pac', data_pac', s_rate,'freqs', phaserange ,'alpha', alpha ,'methodpac', method, 'nfreqs1',nfreqsphase, 'nfreqs2',nfreqsamp,'freqs2',amprange,'winsize',s_rate, 'nboot', nsurrogate);
    elapsedtime = toc;
    fprintf('Elapsed Time for Mean Vector Length Method method %.f s \n', elapsedtime);
    
    % Show plots
    if n_trials == 1,
        eeg_visualize_pac(mvlmipacstruct, 'phasefreq', fm, 'ampfreq', fc, 'abspacval', 1, 'plotall', 1);
    else
        eeg_visualize_pac(mvlmipacstruct, 'time', timepointpac, 'phasefreq', fm, 'ampfreq', fc, 'abs_pacval', 1, 'plot_all', 1, 'nslices', 10);
        eeg_visualize_pac(mvlmipacstruct, 'time', timepointnopac, 'phasefreq', fm, 'ampfreq', fc, 'plot_mvlmi_dist', 1, 'plot_mvlmi_comp', 1, 'plot_comod_time', 1, 'plot_surr', 1, 'abs_pacval', 1);
    end   
end

%% Compute PAC KL
if compute_kl,
    method = 'klmi';
    tic;
    [klmipacval, ~, ~, ~, ~, ~,~,~, klmipacstruct] = eeg_pac(data_pac', data_pac', s_rate,'freqs', phaserange ,'alpha', alpha ,'methodpac', method, 'nfreqs1',nfreqsphase, 'nfreqs2',nfreqsamp,'freqs2',amprange,'winsize',s_rate, 'nbinskl', nbinskl, 'nboot', nsurrogate);
    elapsedtime = toc;
    fprintf('Elapsed Time for Kullback Leibler method %.f s \n', elapsedtime);
    
    % Show plots
    if n_trials == 1,
        eeg_visualize_pac(klmipacstruct, 'phasefreq', fm, 'ampfreq', fc, 'abspacval', 1, 'plotall', 1);
    else
        eeg_visualize_pac(klmipacstruct, 'time', timepointpac, 'phasefreq', fm, 'ampfreq', fc, 'abs_pacval', 1, 'plot_all', 1);
        eeg_visualize_pac(klmipacstruct, 'time', timepointnopac, 'phasefreq', fm, 'ampfreq', fc, 'plot_kl', 1, 'plot_comod_time', 1, 'plot_surr', 1);
    end
end

%% Compute PAC GLM
if compute_glm,
    method = 'glm';
    tic;
    [glmpacval, ~, ~, ~, ~, ~,~,~, glmpacstruct] = eeg_pac(data_pac', data_pac', s_rate,'freqs', phaserange ,'alpha', alpha ,'methodpac', method, 'nfreqs1',nfreqsphase, 'nfreqs2',nfreqsamp,'freqs2',amprange,'winsize',s_rate, 'nboot', nsurrogate);
    elapsedtime = toc;
    fprintf('Elapsed Time for GLM method %.f s \n', elapsedtime);
    
    % Show plots
    if n_trials == 1,
        eeg_visualize_pac(glmpacstruct, 'phasefreq', fm, 'ampfreq', fc, 'abs_pacval', 1, 'plot_all', 1);
    else
        eeg_visualize_pac(glmpacstruct, 'time', timepointpac, 'phasefreq', fm, 'ampfreq', fc, 'abs_pacval', 1, 'plot_all', 1);
        eeg_visualize_pac(glmpacstruct, 'time', timepointnopac, 'phasefreq', fm, 'ampfreq', fc, 'plot_comod_time', 1, 'plot_surr', 1);
    end
end

