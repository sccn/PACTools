% pop_pacprecomp() - precompute PAC measure for a given set of channels or
%                    channels or IC clusters
% Usage:    
%                >> [STUDY, ALLEEG] = pop_pacprecomp(STUDY, ALLEEG); % pop up interactive window
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%
% Outputs:
%   STUDY        - the input STUDY set with added pre-clustering data for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding 
%                  pre-clustering data (pointers to .mat files that hold cluster measure information).
%
% Authors: Ramon Martinez-Cancino,  SCCN, UCSD 2020
%          Arnaud Delorme, SCCN, UCSD 2020
%
% See also: std_pacprecomp()

% Copyright (C) 2020  Ramon Martinez-Cancino, INC, SCCN
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

function [STUDY, ALLEEG, com] = pop_pacprecomp(STUDY, ALLEEG, varargin)
com = '';

% Checking inputs
if nargin < 2
    error('pop_pacprecomp(): Both ALLEEG and STUDY structures are required');
end

if isempty(ALLEEG)
    error('STUDY contains no datasets');
end

%%
% Defaults for GUI
freqs1         = [4 15];
freqs2         = [ALLEEG(1).srate/4-20 ALLEEG(1).srate/4-2];
nfreqs1        = length(freqs1(1):2:freqs1(2));
nfreqs2        = length(freqs2(1):3:freqs2(2)) ;
method_list = {'mvlmi','klmi', 'glm', 'plv', 'instmipac', 'ermipac'} ;

[g,pacargs] = finputcheck(varargin, ...
    { 'datatype'      'string'       {'chan', 'comp'}        'chan';
    'freqs1'          'real'         []                      [4 12]; % Remove default values
    'freqs2'          'real'         []                      [20 60];  % Remove default values
    'freqrange1'      'real'         []                      freqs1;
    'freqrange2'      'real'         []                      freqs2;
    'nfreqs1'         'real'         []                      nfreqs1;
    'nfreqs2'         'real'         []                      nfreqs2;
    'nfreqs2'         'real'         []                      nfreqs2;
    'dataindx'        'integer'      []                      1;
    'method'          'string'       method_list            'glm';
    'freqscale'       'string'       { 'linear','log' }     'linear';
    'compflag'        'string'       {'local', 'nsg'}       'local';
    'nsgopt'          'cell'         {}                      {}}, 'pop_pacprecomp','ignore');

if isempty(g.freqs1), freqs1 = compfreq(g.freqrange1, g.nfreqs1, g.freqscale); end
if isempty(g.freqs2), freqs2 = compfreq(g.freqrange2, g.nfreqs2, g.freqscale); end

if nargin < 3
    pacargs = '';
    % Checking if NSGPORTAL Toolbox is set
    nsginstalled_flag = 1;
    try
        nsg_info;  % From here, get information on where to create the temporary file
    catch
        nsginstalled_flag = 0;
    end

    % Checking if clusters are present
    if length(STUDY.cluster) == 1
        datatype_list = {'Channels'} ;
    else
        datatype_list = {'Channels','Component Clusters'} ;
    end
    
    % Datatype for GUI
    if strcmp(g.datatype, 'chan')
        datatypeindx = 1;
    else
        datatypeindx = 2;
    end
    
    % Method for GUI
    method_indx = find(~cellfun(@isempty,strfind(method_list, g.method)));
    
    method_listgui = {'Mean vector length modulation index (Canolty et al.)',...
                      'Kullback-Leibler modulation index (Tort et al.)',...
                      'General linear model (Penny et al.)',...
                      'Phase Locking Value (Lachaux et al.)',...
                      'Instantaneous MIPAC (Martinez-Cancino et al.)',...
                      'Event related MIPAC (Martinez-Cancino et al.)'};
    guititle = 'Estimate cross-frequency coupling for the STUDY -- pop_pacprecomp()';
    callback_chkcbx_logphs = 'set(findobj(''tag'',''chckbx_logamp''), ''value'', get(findobj(''tag'',''chckbx_logphs''),''value''))';
    callback_chkcbx_logamp = 'set(findobj(''tag'',''chckbx_logphs''), ''value'', get(findobj(''tag'',''chckbx_logamp''),''value''))';
    
    cbdataindx =   ['if get(findobj(''tag'',''datatype''),''value'') == 1,' ...
                    '        [datindx,datstr] = pop_chansel({ALLEEG(1).chanlocs.labels});' ...
                    '   else,' ...
                    '         [datindx,datstr] =  pop_chansel({STUDY.cluster.name});' ...
                    'end;' ...
                    'if ~isempty(datstr)' ...
                    '      set(findobj(''tag'', ''edit_chanind''), ''string'', datstr);' ...
                    'end;'...
                    'set(findobj(''tag'', ''pop_pacprecomp_gui''), ''userdata'', datindx);' ];
                cbpopupdata = 'set(findobj(''tag'', ''edit_chanind''), ''string'', '''');';
    
     nsgcheck = ['if ' num2str(nsginstalled_flag) ',if get(findobj(''tag'',''chckbx_nsgt''),''value''),' ...
                        'set(findobj(''tag'',''nsgopt''),''enable'', ''on'');'...
                     'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off'');end;'...
                   'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off''); set(findobj(''tag'',''chckbx_nsgt''),''value'',0); end'];
    
    geometry =  { {2.4 9 [0    0]  [1    1]}  {2.4 9 [0.6   0] [0.75  1]}...
                  {2.4 9 [0    1]  [1    1]}  {2.4 9 [0.6   1] [1.57  1]} {2.4 9 [2.1 1]   [0.3 1] }...
                  {2.4 9 [0.6  2]  [0.58 1]}  {2.4 9 [1.18 2]  [0.62 1]}  {2.4 9 [1.8 2] [0.45 1] } ...
                  {2.4 9 [0.18 3]  [1    1]}  {2.4 9 [0.6   3] [0.58  1]} {2.4 9 [1.18 3]  [0.62 1]}  {2.4 9 [1.87 3] [0.45 1]} ...
                  {2.4 9 [0.18 4]  [1    1]}  {2.4 9 [0.6   4] [0.58  1]} {2.4 9 [1.18 4]  [0.62 1]}  {2.4 9 [1.87 4] [0.45 1]} ...
                  {2.4 9 [0    5]  [1    1]}  {2.4 9 [0.6   5] [1.8  1]}...
                  {2.4 9 [0    6]  [1    1]}  {2.4 9 [0.6   6] [1.8  1]}...
                  {2.4 9 [0    7]  [1    1]}  {2.4 9 [0.6   7] [1 1]}...
                  {2.4 9 [0    8]  [1    1]}  {2.4 9 [0.6   8] [1.8 1]}};
    
    
    userdata = [];
    uilist = {{'style' 'text' 'string' 'Data type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' datatype_list 'tag' 'datatype' 'value' datatypeindx 'callback' cbpopupdata}...
              {'style' 'text' 'string' 'Clusters/Chan indices' 'fontweight' 'bold' } {'style' 'edit' 'string' ' ' 'tag' 'edit_chanind'}  { 'style' 'pushbutton' 'string' '...' 'callback' cbdataindx } ...
              {'style' 'text' 'string' 'Freq range [lo hi] (Hz)' 'fontweight' 'normal'} {'style' 'text' 'string' '# Frequencies' 'fontweight' 'normal'} {'style' 'text' 'string' 'Log-scaling' 'fontweight' 'normal'} ...
              {'style' 'text' 'string' 'Phase data' 'fontweight' 'normal' 'tag' 'data1'} ...
              {'style' 'edit' 'string' num2str(g.freqrange1)         'tag' 'freq1'}...
              {'style' 'edit' 'string' num2str(g.nfreqs1)        'tag' 'nfreqs1'}...
              {'style' 'checkbox' 'tag' 'chckbx_logphs' 'callback' callback_chkcbx_logphs 'value' 1}...
              {'style' 'text' 'string' 'Amp data ' 'fontweight' 'normal' 'tag' 'data2'} ...
              {'style' 'edit' 'string' num2str(g.freqrange2)        'tag' 'freq2'}...
              {'style' 'edit' 'string' num2str(g.nfreqs2)        'tag' 'nfreqs2'}...
              {'style' 'checkbox' 'tag' 'chckbx_logamp' 'callback' callback_chkcbx_logamp 'value' 1}...
              {'style' 'text' 'string' 'PAC method' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' method_listgui 'value' method_indx 'tag' 'method' }...
              {'style' 'text' 'string' 'Optional inputs' 'fontweight' 'bold'} {'style' 'edit' 'string' ' ' 'tag' 'edit_optinput'}...
              {'style' 'text' 'string' 'Compute on NSG' 'fontweight' 'bold'} {'style' 'checkbox' 'tag' 'chckbx_nsgt' 'callback' nsgcheck 'value' fastif(strcmpi(g.compflag,'nsg'), 1, 0)}...
              {'style' 'text' 'string' 'NSG options'}    {'style' 'edit'       'string' ' '           'tag' 'nsgopt' 'enable' fastif( strcmpi(g.compflag,'nsg'), 'on','off')}};
    
    [out_param dataindx tmp res] = inputgui('title', guititle, 'geom', geometry, 'uilist',uilist,'eval', 'set(gcf,''tag'', ''pop_pacprecomp_gui'')', 'helpcom','pophelp(''pop_pacprecomp'');');
    if isempty(out_param), return; end
    
     % Update g structure
     g.datatype = res.datatype;
     if ~isempty(res.freq1), g.freqrange1   = str2num(res.freq1); else, return; end
     if ~isempty(res.freq2), g.freqrange2   = str2num(res.freq2); else, return; end
     if res.chckbx_logphs,   g.freqscale    = 'log';              else,  g.freqscale = 'linear'; end
     g.dataindx = {ALLEEG(1).chanlocs(dataindx).labels};
     g.method = method_list{res.method};
     g.freqs1 = compfreq(g.freqrange1, g.nfreqs1, g.freqscale);
     g.freqs2 = compfreq(g.freqrange2, g.nfreqs2, g.freqscale);   
     
     % Update fields present here
     tmpparams = eval( [ '{' res.edit_optinput '}' ] );
     for i =1:2:length(tmpparams)
             g.(tmpparams{i}) = tmpparams{i+1};
     end
     
     if res.chckbx_nsgt
         g.nsgopt = eval(['{''nsgflag'', 1,' res.nsgopt '}']);
     else
         g.nsgopt = {'nsgflag', 0};
     end  
end

 % Removing fields not needed in std_pacprecomp
fields2rm = {'freqrange1', 'freqrange2', 'nfreqs1', 'nfreqs2', 'freqscale'};
g = rmfield(g,fields2rm);
opt = struct2args(g);

% Calling std_precomp
[STUDY, ALLEEG] = std_pacprecomp(STUDY, ALLEEG, freqs1, freqs2, g.datatype, g.dataindx, opt{:});

% Command call
com = sprintf('[STUDY, ALLEEG] = std_pacprecomp(STUDY, ALLEEG, %s);', vararg2str(opt(3:end)));


function freqs = compfreq(freqrange, nfreqs, freqscale)
if strcmpi(freqscale, 'log')
        freqs = linspace(log(freqrange(1)), log(freqrange(end)), nfreqs);
        freqs = exp(freqs);
    else
        freqs = linspace(freqrange(1), freqrange(2), nfreqs);
    end