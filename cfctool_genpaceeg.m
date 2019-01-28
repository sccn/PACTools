%% Generating single trial PAC signal
% Signal parameters
fc            = 60;   % Frequency of the carrier wave. Recommended: [60]
fm            = 8;    % Frequency of the modulating wave. Recommended: [60]
max_time      = 7;    % Time of the simulation. Recommended: single trial [6], multiple trials [2]
n_trials      = 300;    % Number of trials. Recommended: single trial [1], multiple trials [1000]
s_rate        = 500;  % Sampling rate. Recommended: single trial [500], multiple trials [250]
padtime       = 0;    % Padding time at the start and end of the signal. 
                      % Should increase with sampling size.
                      % Recommended: 0.5
snrval        = 10;    % Signal to noise ratio. Recommended: [5]
maxshift      = 10;   % Maximum shifts (jitter) injected into the data.Recommended: single trial [1], multiple trials [10]
nsegm         = 5;    % Number of segments in the data. Each segment is a block with our without modulation

% Do not modify the code beyond this point
%% Generate 1D data
tlimits      = [0 max_time];
plot_flag    = 1;
[data_singletrial,t_outtmp]  = cfctool_simpac(fc,fm,tlimits,'Ac',5,'Am',1,'fs',s_rate,'cpfunc','block','blockamp',1, 'nsegm',nsegm,'plot_flag', plot_flag,'padtime',padtime,'snr',snrval,'m',1);
                                                                                 
%%
EEG = eeg_emptyset;
EEG.setname  = 'simPACdata';
EEG.filename = 'strial_sim_pac_ph8amp60';
EEG.nbchan   = 1;
EEG.trials   = 1;
EEG.pnts     = length(data_singletrial);
EEG.xmin     = t_outtmp(1);
EEG.xmax     = t_outtmp(end);
EEG.data     = data_singletrial;
EEG.srate    = s_rate;
EEG.times    = t_outtmp;
EEG = pop_saveset(EEG, 'filename','strial_sim_pac_ph8amp60' ,'filepath',pwd,'savemode', 'onefile' );

%% Generate 1D EEG data (single trial, 2 channel)
tlimits      = [0 max_time];
plot_flag    = 1;
[data_singletrial,t_outtmp]  = cfctool_simpac(fc,fm,tlimits,'Ac',5,'Am',1,'fs',s_rate,'cpfunc','block','blockamp',1, 'nsegm',nsegm,'plot_flag', plot_flag,'padtime',padtime,'snr',snrval,'m',1);
                                                                                 
%%
EEG = eeg_emptyset;
EEG.setname  = 'simPACdata';
EEG.filename = 'strial_sim_pac_ph8amp60';
EEG.nbchan   = 2;
EEG.trials   = 1;
EEG.pnts     = length(data_singletrial);
EEG.xmin     = t_outtmp(1);
EEG.xmax     = t_outtmp(end);
EEG.data     = [data_singletrial;awgn(data_singletrial,10)] ;
EEG.srate    = s_rate;
EEG.times    = t_outtmp;
EEG = pop_saveset(EEG, 'filename','strial_2chan_sim_pac_ph8amp60' ,'filepath',pwd,'savemode', 'onefile' );


%% Generate 2d EEGdata (multiple trials, 1 channel)

max_time      = 2;            % Time of the simulation. Recommended: single trial [6], multiple trials [2]
nfreqsphase   = 1;
nfreqsamp     = 1;

tlimits = [-1 max_time];
s_rate  = 50; 
snrval  = 10;     
s = randi(maxshift, n_trials,1);

for i=1:n_trials
     plot_flag = i == 1;
     [datatmp,t_outtmp]  = cfctool_simpac(fc,fm,tlimits,'Ac',5,'Am',1,'fs',s_rate,'cpfunc','block','blockamp',1, 'nsegm',nsegm,'plot_flag', plot_flag,'padtime',padtime,'snr',snrval,'m',1);
     datatmp = [datatmp(s(i):end) datatmp(1:s(i)-1)];
     data_pac(1,:,i) = [datatmp(end-ceil(maxshift/2):end) datatmp(1:end-ceil(maxshift/2)-1)]; 
end

%
EEG = eeg_emptyset;
EEG.setname  = 'simPACdata_multtrial';
EEG.filename = 'mult_trial_sim_pac_ph8amp60';
EEG.nbchan   = 1;
EEG.trials   = n_trials;
EEG.pnts     = length(datatmp);
EEG.xmin     = t_outtmp(1);
EEG.xmax     = t_outtmp(end);
EEG.data     = data_pac ;
EEG.srate    = s_rate;
EEG.times    = t_outtmp;
% EEG = pop_saveset(EEG, 'filename','mult_trial_2chan_sim_pac_ph8amp60' ,'filepath',pwd,'savemode', 'onefile' );


