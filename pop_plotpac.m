% Inputs:
%    EEG      - EEG structure or pac structure
%
%   Optional inputs
%       'time'              - Time point(s) to plot in seconds. Only one  
%                             time point is used for all the plots except
%                             plotcomodt. Default [1]                
%       'phasefreq'         - Frequency of the phase to plot. The plots will
%                             use the closest frequency available and the
%                             median frequency by default                           
%       'ampfreq'           - Frequency of the amplitude to plot. The plots will
%                             use the closest frequency available and the
%                             median frequency be default.
%       'plotcomod'         - Logical. Plot the comodulogram. 
%                             Default [1] if no other plot flags given,
%                             otherwise [0]
%       'plotcomodt'        - Logical. Plot the comodulogram at the given 
%                             time point(s). Several time points can be used.
%                             Default [0]
%       'plotampt'          - Logical. Plot the modulation index in time
%                             and in amplitudes at a given phase frequency.
%                             Default [0]
%       'plotphaset'        - Logical. Plot the modulation index in time
%                             and in phases at a given amplitude frequency.
%                             Default [0]
%       'plotsurrdist'      - Logical. Plot the surrogate distribution for 
%                             a given time point, amplitude frequency, 
%                             and phase frequency. Default [0]                          
%       'plotkl'            - Logical. Plot Kullback-Leibler modulation 
%                             index plot phase amplitude for a given time point,
%                             amplitude frequency, and phase frequency.                          
%                             Default [0]   
%       'plotmvl'           - Logical. Plot the distribution of the composites 
%                             of the mean vector length modulation index method
%                             for a given time point, amplitude frequency,
%                             and phase frequency. Default [0]
%       'arrowweight'       - Scalar used to multiply the ploted mean 
%                             vector length arrow. Default [0]
%       'nbinsmvl'          - Number of bins per phase bin and half the number 
%                             of phase bins to use in the composite distribution 
%                             of the mean vector length modulation index plot. 
%                             Default [36]
%       'normcomposite'     - Logical. Normalize the amplitude of the composites 
%                             values of the mean vector length modulation index. 
%                             Default [0].
%       'abspacval'         - Use the absolute value of the phase amplitude
%                             coupling value. Default [0]
%       'plotall'           - Logical. Plot all of the possible plots with
%                             the data and given parameters. Default [0]
%       'alphadata'         - Transparency to be used for to mask of
%                             significance in the comodulogram plot.
%                             Default [0.3]
%       'phasechanindx'     - When EEG structure is provided, will take
%                             this index to retreive and reconstruct the
%                             'pacstruct'
%       'ampchanindx'       - When EEG structure is provided, will take
%                             this index to retreive and reconstruct the
%                            'pacstruct'
function  hfig  = pop_plotpac(EEG, varargin)
hfig = [];
% Checking EEG structure
if nargin == 0
    help pop_plotpac;
    return;
end
if isempty(EEG)
    fprintf(2,'pop_plotpac error : Empty input provided \n');
    return;
end
EEG = eeg_checkset(EEG);

 % Check pac field on EEG
if ~isfield(EEG,'etc') || ~isfield(EEG.etc,'eegpac') || isempty(EEG.etc.eegpac)
    fprintf(2,'pop_plotpac error : PAC must be computed first \n');
    return;
end

%% Inputs
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('pop_plotpac() error: calling convention {''key'', value, ... } error'); return;
end;

timeflag  = 0;
phaseflag = 0;
ampflag   = 0;

try g.time;                timeflag    = 1; catch, g.time           = 1;              end;
try g.phasefreq;           phaseflag   = 1; catch, g.phasefreq      = [];               end;  % Default set once pactruct is defined
try g.ampfreq;             ampflag     = 1; catch, g.ampflag        = [];                                  end;  % Default set once pactruct is defined
try g.plotcomod;           catch, g.plotcomod      = 0;              end;
try g.plotcomodt;          catch, g.plotcomodt     = 0;              end;
try g.plotampt;            catch, g.plotampt       = 0;              end;
try g.plotphaset;          catch, g.plotphaset     = 0;              end;
try g.plotsurrdist;        catch, g.plotsurrdist   = 0;              end;
try g.plotkl;              catch, g.plotkl         = 0;              end;
try g.plotmvl;             catch, g.plotmvl        = 0;              end;
try g.normcomposite;       catch, g.normcomposite  = 0;              end;
try g.nbinsmvl;            catch, g.nbinsmvl       = 36;             end;
try g.abspacval;           catch, g.abspacval      = 0;              end;
try g.plotall;             catch, g.plotall        = 0;              end;
try g.alphadata;           catch, g.alphadata      = 0.3;            end;
try g.arrowweight;         catch, g.arrowweight    = 1;              end;

try g.phasechanindx;       catch, g.phasechanindx  = 1;             end;
try g.ampchanindx;         catch, g.ampchanindx    = 1;             end;

% Set default phase and/or amplitude frequency
if ~phaseflag
    g.phasefreq = median(EEG.etc.eegpac.freqs_phase);
end
if ~ampflag
    g.ampfreq = median(EEG.etc.eegpac.freqs_amp);
end

% Determine if the data is single strial or multiple triale
if numel(size(EEG.etc.eegpac.datapac{1}.pacval)) == 2,
    s_trial = true;
else
    s_trial = false;
end

if isempty(varargin);
    %--
    g.plotcomod = 1;
    g.plotampt  = 1;
    g.plotampt   = 1;
    g.plotphaset = 1;
    g.plotcomodt = 1;
    plotampt_enable   = 'on'; plotphaset_enable = 'on'; plotcomodt_enable = 'on'; 
    if s_trial,
        g.plotampt   = 0; plotampt_enable   = 'off';
        g.plotphaset = 0; plotphaset_enable = 'off';
        g.plotcomodt = 0; plotcomodt_enable = 'off';
    end
    if ~isempty(EEG.etc.eegpac.datapac{g.phasechanindx,g.ampchanindx}.surrogate_pac),
        g.plotsurrdist = 1;
    end
    switch EEG.etc.eegpac.method
        case 'mvlmi'
            g.plotmvl = 1;
        case 'klmi'
            g.plotkl = 1;
    end
    %--
    phaseindxval = 1;
    ampindxval   = 1;
    checkboxabler_callback = '';
    plot_callback = ['hfig = gcbf; hfig = pop_plotpac(EEG,  ''time'',             str2num(get(findobj(hfig,''tag'',''edit_comodtime''),''string'')),'...
                                                            '''phasefreq'',       str2num(get(findobj(hfig,''tag'',''edit_freqphase''),''string'')),'...
                                                            '''ampfreq'',         str2num(get(findobj(hfig,''tag'',''edit_freqamp''),  ''string'')),'...
                                                            '''plotcomod'',       get(findobj(hfig,''tag'',''chckbox_comod''),''value''),'...
                                                            '''plotcomodt'',      get(findobj(hfig,''tag'',''chckbox_comodt''),''value''),'...
                                                            '''plotampt'',        get(findobj(hfig,''tag'',''chckbox_phvsall''),''value''),'...
                                                            '''plotphaset'',      get(findobj(hfig,''tag'',''chckbox_ampvsall''),''value''),'...
                                                            '''plotsurrdist'',    get(findobj(hfig,''tag'',''chckbox_surrdist''),''value''),'...
                                                            '''plotkl'',          get(findobj(hfig,''tag'',''chckbox_klhist''),''value''),'...
                                                            '''plotmvl'',         get(findobj(hfig,''tag'',''chckbox_mvhist''),''value''),'...
                                                            '''normcomposite'','  num2str(g.normcomposite) ','...
                                                            '''nbinsmvl'',        str2num(get(findobj(hfig,''tag'',''edit_mvbins''),''string'')),',...
                                                            '''abspacval'','      num2str(g.abspacval) ','...
                                                            '''alphadata'','      num2str(g.alphadata) ','...
                                                            '''arrowweight'','    num2str(g.arrowweight) ','...
                                                            '''phasechanindx'',   str2num(get(findobj(hfig,''tag'',''edit_phaseindx''),''string'')),'...
                                                            '''ampchanindx'',     str2num(get(findobj(hfig,''tag'',''edit_ampindx''),''string'')))'];
    % Building the GUI
    guititle = 'pop_plotpac()';
    uilist = { ...
        { 'style', 'text', 'string','PAC Component / Channel Indices','fontweight' 'bold' },...
        { 'style', 'text', 'string','Phase Signal'},     {'style' 'edit' 'string' '1' 'callback' checkboxabler_callback 'tag' 'edit_phaseindx'},{ 'style', 'text', 'string', 'Amplitude Signal'},{'style' 'edit' 'string' '1' 'tag' 'ampchanindx' 'callback' checkboxabler_callback 'tag' 'edit_ampindx'}...
        { 'style', 'text', 'string','PAC Plots','fontweight' 'bold'},{ 'style', 'text', 'string','PAC Plots Settings','fontweight' 'bold'}...
        { 'style', 'text', 'string','Comodulogram'},              {'Style', 'checkbox', 'value', g.plotcomod,   'tag' 'chckbox_comod'  },...
        { 'style', 'text', 'string','Temporal Comodulogram'},     {'Style', 'checkbox', 'value', g.plotcomodt,  'tag' 'chckbox_comodt'  'Enable'  plotcomodt_enable }, {'style', 'text', 'string', 'Time (s)'},            {'style' 'edit' 'string' '0.5 : 0.5 : 1.5' 'tag' 'edit_comodtime'} ...
        { 'style', 'text', 'string','Phase Vs all Amplitudes'},   {'Style', 'checkbox', 'value', g.plotampt,    'tag' 'chckbox_phvsall' 'Enable' plotampt_enable },   {'Style', 'text', 'string', 'Phase Freq. Range [low high] (Hz)'},     {'style' 'edit' 'string' '2' 'tag' 'edit_freqphase'} ...
        { 'style', 'text', 'string','Amplitude Vs all Phases'},   {'Style', 'checkbox', 'value', g.plotphaset,  'tag' 'chckbox_ampvsall' 'Enable' plotphaset_enable }, {'Style', 'text', 'string', 'Amplitude Freq. Range [low high] (Hz)'},{'style' 'edit' 'string' '60' 'tag' 'edit_freqamp'} ...
        { 'style', 'text', 'string','Mean-Vector Histogram'},     {'Style', 'checkbox', 'value', g.plotmvl,     'tag' 'chckbox_mvhist'},   {'Style', 'text','string' '# bins'},                {'style' 'edit' 'string' '15' 'tag' 'edit_mvbins'} ...
        { 'style', 'text', 'string','Kullback-Leibler Histogram'},{'Style', 'checkbox', 'value', g.plotkl,      'tag' 'chckbox_klhist'}...
        { 'style', 'text', 'string','Surrogate Distribution'},    {'Style', 'checkbox', 'value', g.plotsurrdist,'tag' 'chckbox_surrdist'} ...
        { 'style', 'text', 'string','Optional eeg_plotpac() arguments(see help)','fontweight' 'bold' },{'style' 'edit' 'string' ' ' 'tag' 'edit_optionalparams'} ...
        { }...
        { 'style', 'pushbutton' , 'string', 'Help', 'callback','pophelp(''pop_plotpac'');' } { 'style', 'pushbutton' , 'string', 'Cancel' 'Callback' 'close(gcbf);'} { 'style', 'pushbutton' , 'string', 'Plot' 'Callback'  plot_callback}...
        };
    
    geometry = {{2.3 14 [0 0]  [1 1]}...
        {2.3 14 [0 1]  [1 1]}     {2.3 14 [0.5 1]  [0.3 1]} ...
        {2.3 14 [0 2]  [1 1]}     {2.3 14 [0.5 2]  [0.3 1]}...
        {2.3 14 [0 3]  [1 1]}     {2.3 14 [1 3]  [1 1]}...
        {2.3 14 [0 4]  [1 1]}     {2.3 14  [0.65 4]  [0.2 1]}...
        {2.3 14 [0 5]  [1 1]}     {2.3 14  [0.65 5]  [0.2 1]} {2.3 14 [1 4] [1 1]}  {2.3 14 [1.8 4] [0.5 1]}...
        {2.3 14 [0 6]  [1 1]}     {2.3 14  [0.65 6]  [0.2 1]} {2.3 14 [1 5] [1 1]}  {2.3 14 [1.8 5] [0.5 1]}...
        {2.3 14 [0 7]  [1 1]}     {2.3 14  [0.65 7]  [0.2 1]} {2.3 14 [1 6] [1 1]}  {2.3 14 [1.8 6] [0.5 1]}...
        {2.3 14 [0 8]  [1 1]}     {2.3 14  [0.65 8]  [0.2 1]} {2.3 14 [1 7] [1 1]}  {2.3 14 [1.8 7] [0.5 1]}...
        {2.3 14 [0 9]  [1 1]}     {2.3 14  [0.65 9]  [0.2 1]}...
        {2.3 14 [0 10]  [1 1]}    {2.3 14  [0.65 10]  [0.2 1]} ...
        {2.3 14 [0 11] [1 1]}     {2.35 14 [1 11] [1.35 1]}...
        {2.3 14 [0 12] [1 1]}...
        {2.3 14 [0 13] [0.5 1]}  {2.3 14 [1.33 13] [0.5 1]}  {2.3 14 [1.8 12] [0.5 1]}};
    
    
    handles = supergui('title', guititle, 'geom', geometry, 'uilist',uilist);
else
    options = {'time',          g.time,...
               'phasefreq',     g.phasefreq,...
               'ampfreq',       g.ampfreq,...
               'plotcomod',     g.plotcomod,...
               'plotcomodt',    g.plotcomodt,...
               'plotampt',      g.plotampt,...
               'plotphaset',    g.plotphaset,...
               'plotsurrdist',  g.plotsurrdist,...
               'plotkl',        g.plotkl,...
               'plotmvl',       g.plotmvl,...
               'normcomposite'  g.normcomposite,...
               'nbinsmvl',      g.nbinsmvl,...
               'abspacval',     g.abspacval,...
               'plotall',       g.plotall,...
               'alphadata',     g.alphadata,...
               'arrowweight',   g.arrowweight,...
               'phasechanindx', g.phasechanindx,...
               'ampchanindx',   g.ampchanindx};
                
    hfig = eeg_visualize_pac(EEG, options{:});
end
end
