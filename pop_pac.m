% pop_pac() - Call GUI to compute cross-frequency-coupling coupling.
%             Second level function to compute CFC by calling eeg_pac.m
% Usage:
%   >>  pac = pop_pac(EEG);
%
% Inputs:
%  EEG          - [Structure] Input dataset as an EEGLAB EEG structure
%  pooldata     - ('channels' || 'component' ). Define if compute CFC in
%                 channel or ICA decomposed data
%  freqs1       - [min max] Range of frequency to consider for the phase values.
%  freqs2       - [min max] Range of frequency to consider for the amplitude values.
%  indexfreqs1  - [Integer]Index of the channel or component to use to in freqs1
%  indexfreqs2  - [Integer]Index of the channel or component to use to in freqs2
%
% Optional inputs:
% Note: Optional parameters of eeg_pac can be provided as input as well
%  method       - {'mvlmi','klmi','glm','instmipac', 'ermipac'}.CFC method .
%                 Default: 'glm'. (currently PAC methods only)
%  nboot        - [Integer] Number of surrogates generated for statistical significance analysis. Default: [200] 
%  alpha        - [Real] Significance threshold. Default: [0.05]
%  bonfcorr     - [1,0] Flag to perform multiple comparisons correction.Default: [0] 
%                 using Bonferroni. Default: [0](do not perform correction)
%  nfreqs1      - [Integer]Number of frequencies in the 'freqs1' range. Default: [1] 
%  nfreqs2      - [Integer]Number of frequencies in the 'freqs2' range. Default: [1] 
%  cleanup      - [1,0] Flag to remove previous results in the EEG
%                 structure. Default: [0] (Do not clean up structure)
%  ptspercent   - [0.01:0.5] Size in percentage of data of the segments to shuffle 
%                 when creating surrogate data. Default: [0.05] 
%  forcecomp    - [0,1] Flag to force (1) or not (0) the computation of PAC
%                 in the case it is detected that the measure has been
%                 computed already or the parameters used are different to the ones saved.
%                 Default:[0]
% Outputs:
%  EEG          - EEGLAB EEG structure. Results of computing CFC are
%                  stored in EEG.etc.pac
%  com          - Command for EEG history
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2019
%
% Copyright (C) 2019  Ramon Martinez-Cancino,INC, SCCN
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

function [EEG,com] = pop_pac(EEG,pooldata,freqs1,freqs2,indexfreqs1,indexfreqs2,varargin)
com = '';

if nargin == 0
    help pop_pac;
    return;
end
if isempty(EEG) || isempty(EEG.data)
    fprintf(2,'pop_pac error : Unsuported input has been provided \n');
    return;
end
EEG = eeg_checkset(EEG);

% get inputs here
[tmp,tmp,dim_trial] = size(EEG.data); clear tmp;

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else
        g = [];
    end
catch
    disp('pop_pac() error: calling convention {''key'', value, ... } error'); return;
end
try g.freq1_trialindx;      catch, g.freq1_trialindx  =  1:dim_trial; end % Opt only for pop_pac
try g.freq2_trialindx;      catch, g.freq2_trialindx  =  1:dim_trial; end % Opt only for pop_pac

try g.nboot;                catch, g.nboot            =  200;         end
try g.alpha;                catch, g.alpha            =  0.05;        end
try g.bonfcorr;             catch, g.bonfcorr         =  0;           end
try g.method;               catch, g.method           =  'glm';       end
try g.nfreqs1;              catch, g.nfreqs1          =  1;           end
try g.nfreqs2;              catch, g.nfreqs2          =  1;           end
try g.cleanup;              catch, g.cleanup          =  0;           end
try g.freqscale;            catch, g.freqscale        =  'log';       end % may be 'linear' or 'log'
try g.ptspercent;           catch, g.ptspercent       =  0.05;        end 
try g.forcecomp;            catch, g.forcecomp        =  0;           end
try g.compflag;             catch, g.compflag         = 'local';      end
try g.runtime;              catch, g.runtime          = 1;            end
try g.jobid;                catch, g.jobid            = ['pacssnsg_' num2str(floor(rand(1)*1000000))]; end

if nargin < 6 
    
    % Closing open GUI and creating a new one
    openfig = findobj('tag', 'pop_pacgui');
    if ~isempty(openfig)
        disp('pop_pac warning: there can be only one pop_pac window, closing old one...')
        close(openfig); 
    end
    
    % Checking NSG
    nsginstalled_flag = 1;
    try
        nsg_info;  % From here, get information on where to create the temporary file
    catch
        nsginstalled_flag = 0;
    end
    
    % Defaults for GUI
    freqs1         = [4 15];
    freqs2         = [EEG.srate/4-20 EEG.srate/4-2];
    nfreqs1        = length(freqs1(1):2:freqs1(2));
    nfreqs2        = length(freqs2(1):3:freqs2(2)) ;
    freq1_dataindx = 1;
    freq2_dataindx = 1;
    
    % ---
    % Checking if ICA
    if isempty(EEG.icawinv)
        datatypel_list = {'Channels'} ;
    else
        datatypel_list = {'Channels','Components'} ;
    end
    
    % Here define other types of CFC
    % Note: In case of adding more methods, add the modality here at the
    % end of the cell array, and then update 'data1_list' and 'data2_list'
    % in the callback 'callback_setdata'
    %%
        cb_popmenu_datatype = 'set(findobj(''tag'', ''freq1_dataindx''), ''string'', '' ''); set(findobj(''tag'', ''freq2_dataindx''), ''string'', '' '');';
    nsgcheck = ['if ' num2str(nsginstalled_flag) ',if get(findobj(''tag'',''chckbx_nsgt''),''value''),' ...
                        'set(findobj(''tag'',''nsgopt''),''enable'', ''on'');'...
                     'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off'');end;'...
                   'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off''); set(findobj(''tag'',''chckbx_nsgt''),''value'',0); end'];
    
%    cfctype_list = {'Phase-Amp'};%, 'Amp.-Amp.', 'Phase-Phase'}; 
    
   cbdataindx1 =   [ 'set(findobj(''tag'', ''freq1_dataindx''), ''string'', '' '');'...
                    'if get(findobj(''tag'',''datatype''),''value'') == 1,' ...
                    '   if isempty(EEG(1).chanlocs),'...
                    '     chanlist = cellfun(@(x) [''Chan'' num2str(x)],num2cell(1:size(EEG.data,1)), ''UniformOutput'',0);'...
                    '     [datindx1,datstr] = pop_chansel(chanlist);'...
                    '   else;'...
                    '        [datindx1,datstr] = pop_chansel({EEG(1).chanlocs.labels});' ...
                    '   end;'...
                    '   else,' ...
                    '         iclist = cellfun(@(x) [''IC'' num2str(x)],num2cell(1:size(EEG.icawinv,2)), ''UniformOutput'',0);'...
                    '         [datindx1,datstr] =  pop_chansel(iclist);' ...
                    'end;' ...
                    'if ~isempty(datstr)' ...
                    '      set(findobj(''tag'', ''freq1_dataindx''), ''string'', datstr);' ...
                    'end;'...
                    'myappdata = getappdata(findobj(''tag'', ''pop_pacgui''),''userdata'');'...
                    'myappdata.freqindx1 = datindx1;'...
                    'setappdata(findobj(''tag'', ''pop_pacgui''), ''userdata'', myappdata);'...
                    'set(findobj(''tag'', ''pop_pacgui''), ''userdata'', myappdata);' ];
                
   cbdataindx2 =   [ 'set(findobj(''tag'', ''freq2_dataindx''), ''string'', '' '');'...
                    'if get(findobj(''tag'',''datatype''),''value'') == 1,' ...
                    '   if isempty(EEG(1).chanlocs),'...
                    '     chanlist = cellfun(@(x) [''Chan'' num2str(x)],num2cell(1:size(EEG.data,1)), ''UniformOutput'',0);'...
                    '     [datindx2,datstr] = pop_chansel(chanlist);'...
                    '   else;'...
                    '        [datindx2,datstr] = pop_chansel({EEG(1).chanlocs.labels});' ...
                    'end;'...
                    '   else,' ...
                    '         iclist = cellfun(@(x) [''IC'' num2str(x)],num2cell(1:size(EEG.icawinv,2)), ''UniformOutput'',0);'...
                    '         [datindx2,datstr] =  pop_chansel(iclist);' ...
                    'end;' ...
                    'if ~isempty(datstr)' ...
                    '      set(findobj(''tag'', ''freq2_dataindx''), ''string'', datstr);' ...
                    'end;'...
                    'myappdata = getappdata(findobj(''tag'', ''pop_pacgui''),''userdata'');'...
                    'myappdata.freqindx2 = datindx2;'...
                    'setappdata(findobj(''tag'', ''pop_pacgui''), ''userdata'', myappdata);'...
                    'set(findobj(''tag'', ''pop_pacgui''), ''userdata'', myappdata);' ];
                
   
    method_listgui = {'Mean vector length modulation index (Canolty et al., 2006)',...
                      'Kullback-Leibler modulation index (Tort et al., 2010)',...
                      'General linear model (Penny et al., 2008)',...
                      'Phase Locking Value (Lachaux et al. 1999)',...
                      'Instantaneous MIPAC (Martinez-Cancino et al., 2019)',...
                      'Event related MIPAC (Martinez-Cancino et al., 2019)'};
                  
    method_list    = {'mvlmi','klmi','glm','plv','instmipac', 'ermipac'};
    guititle = 'Estimate event-related cross-frequency coupling -- pop_pac()'; 
    callback_chkcbxstat = ['label_statsate = {''(on)'', ''(off)''};'... 
                           'set(findobj(''tag'', ''label_statstate''), ''string'', label_statsate{fastif(get(findobj(gcf,''tag'', ''chckbx_stat''), ''value''),1,2)});'...
                           'if get(findobj(''tag'',''chckbx_stat''),''value''),' ...
                           'set(findobj(''tag'',''nsurrogates''),''enable'', ''on'');'...
                           'set(findobj(''tag'',''nblocks_edit''),''enable'', ''on'');'...
                           'set(findobj(''tag'',''pvalue''),''enable'', ''on'');'...
                           'set(findobj(''tag'',''bonfcorr''),''enable'', ''on'');'...
                           'else, set(findobj(''tag'',''nsurrogates''),''enable'', ''off'');'...
                           'set(findobj(''tag'',''nblocks_edit''),''enable'', ''off'');'...
                           'set(findobj(''tag'',''pvalue''),''enable'', ''off'');'...
                           'set(findobj(''tag'',''bonfcorr''),''enable'', ''off'');end;'];
                     
    callback_chkcbx_logphs = 'set(findobj(''tag'',''chckbx_logamp''), ''value'', get(findobj(''tag'',''chckbx_logphs''),''value''))';
    callback_chkcbx_logamp = 'set(findobj(''tag'',''chckbx_logphs''), ''value'', get(findobj(''tag'',''chckbx_logamp''),''value''))';
    
    guiheight = 12;
    guiwidth = 2.6;
    if nsginstalled_flag
        nsgdefaultopt = ['''runtime'',' num2str(g.runtime) ',''jobid'','''  g.jobid ''''];
        nsgmenugeom = {{guiwidth guiheight [0 10]    [1 1]}  {guiwidth guiheight [0.5 10] [0.5  1]}...
            {guiwidth guiheight [0.1 11]    [1 1]}  {guiwidth guiheight [0.5 11] [2.1 1]}...
            {guiwidth guiheight [0 12]    [1 1]}};
        nsgmenuuilist = {{'style' 'text' 'string' 'Compute via NSG' 'fontweight' 'bold'} {'style' 'checkbox' 'tag' 'chckbx_nsgt' 'callback' nsgcheck 'value' fastif(strcmpi(g.compflag,'nsg'), 1, 0)}...
            {'style' 'text' 'string' 'NSG options'}    {'style' 'edit'       'string'  nsgdefaultopt           'tag' 'nsgopt' 'enable' fastif( strcmpi(g.compflag,'nsg'), 'on','off')}...
            {}};
    else
        guiheight = 9;
        nsgmenugeom =  {{guiwidth guiheight [0 10]    [1 1]}};
        nsgmenuuilist = {{}};
    end
    
    userdata = [];
    geometry = { {guiwidth guiheight [0 0]     [1 1]}  {guiwidth guiheight [0.3  0] [0.8 1] }...%  {guiwidth 9 [1.58 0]  [1 1]}    {guiwidth 9 [1.9 0] [0.83 1] }...
                                                  {guiwidth guiheight [0.6  1] [0.58 1]}                                    {guiwidth guiheight [1.37 1] [0.55 1]} {guiwidth guiheight [1.9 1] [0.38 1]} {guiwidth guiheight [2.25 1] [0.4 1]}...
                 {guiwidth guiheight [0.1 2]  [1 1]}  {guiwidth guiheight [0.6  2] [0.6 1]} {guiwidth guiheight [1.12  2] [0.24 1]}  {guiwidth guiheight [1.37 2] [0.55 1]} {guiwidth guiheight [1.9 2] [0.38 1]} {guiwidth guiheight [2.32 2] [0.3 1]}...
                 {guiwidth guiheight [0.1 3]  [1 1]}  {guiwidth guiheight [0.6  3] [0.6 1]} {guiwidth guiheight [1.12  3] [0.24 1]}  {guiwidth guiheight [1.37 3] [0.55 1]} {guiwidth guiheight [1.9 3] [0.38 1]} {guiwidth guiheight [2.32 3] [0.3 1]}...
                 {guiwidth guiheight [0 4]     [1 1]}  {guiwidth guiheight [0.58  4] [2.03 1]}...
                 {guiwidth guiheight [0 5]     [1 1]}  {guiwidth guiheight [0.6  5]  [2  1]}...
                 {guiwidth guiheight [0 6]     [1 1]}  {guiwidth guiheight [0.6  6] [2.1  1]}...
                 {guiwidth guiheight [0.15 7]  [1 1]}  {guiwidth guiheight [1 7] [0.5  1]}  {guiwidth guiheight [1.75  7] [1 1]} {guiwidth guiheight [2.1  7] [0.5 1]}...
                 {guiwidth guiheight [0.15 8]  [1 1]}  {guiwidth guiheight [1 8] [0.5  1]}...
                 {guiwidth guiheight [0.15 9]  [1 1]}  {guiwidth guiheight [1 9] [0.5  1]}...
                 nsgmenugeom{:}};
             
    uilist = {{'style' 'text' 'string' 'Data type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' datatypel_list 'tag' 'datatype' 'value' 1 'callback' cb_popmenu_datatype}...% {'style' 'text' 'string' 'CFC type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' cfctype_list 'tag' 'datatyppac' 'value' 1 'callback' callback_setdata}...% 
              {'style' 'text' 'string' 'Select ICs|Channels' 'fontweight' 'normal'} {'style' 'text' 'string' 'Freq range [lo hi] (Hz)' 'fontweight' 'normal'} {'style' 'text' 'string' '# Frequencies' 'fontweight' 'normal'} {'style' 'text' 'string' 'Log-scale' 'fontweight' 'normal'} ...
              {'style' 'text' 'string' 'Phase data (Low-freq)' 'fontweight' 'normal' 'tag' 'data1'} {'style' 'edit' 'string' ' ' 'tag' 'freq1_dataindx'} { 'style' 'pushbutton' 'string' '...' 'callback' cbdataindx1 } ...
                                                                                                {'style' 'edit' 'string' num2str(freqs1)         'tag' 'freq1'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs1)        'tag' 'nfreqs1'}...
                                                                                                {'style' 'checkbox' 'tag' 'chckbx_logphs' 'callback' callback_chkcbx_logphs 'value' 1}...
              {'style' 'text' 'string' 'Amplitude data (Hi-freq)' 'fontweight' 'normal' 'tag' 'data2'} ...
                                                                                                {'style' 'edit' 'string' ' ' 'tag' 'freq2_dataindx'} { 'style' 'pushbutton' 'string' '...' 'callback' cbdataindx2 } ......
                                                                                                {'style' 'edit' 'string' num2str(freqs2)        'tag' 'freq2'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs2)        'tag' 'nfreqs2'}...
                                                                                                {'style' 'checkbox' 'tag' 'chckbx_logamp' 'callback' callback_chkcbx_logamp 'value' 1}...
              {'style' 'text' 'string' 'PAC method to compute' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' method_listgui 'tag' 'method'}...
               {'style' 'text' 'string' 'Optional pop_pac inputs' 'fontweight' 'bold'} {'style' 'edit' 'string' ' ' 'tag' 'edit_optinput'}...
              {'style' 'text' 'string' 'Significance testing' 'fontweight' 'bold'} {'style' 'checkbox' 'tag' 'chckbx_stat' 'callback' callback_chkcbxstat 'value' 0}...
              {'style' 'text' 'string' '# surrogates to compute' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.nboot) 'tag' 'nsurrogates' 'enable' 'off'}  {'style' 'text' 'string' '# data blocks' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(1/g.ptspercent) 'tag' 'nblocks_edit' 'enable' 'off'} ...
              {'style' 'text' 'string' 'Significance threshold (0<p<1)' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.alpha) 'tag' 'pvalue' 'enable' 'off'}... 
              {'style' 'text' 'string' 'Correct for multiple comparisons (BNF)' 'callback' 'close(gcbf);' } {'style' 'checkbox' 'tag' 'bonfcorr' 'enable' 'off'}...
              nsgmenuuilist{:}};
          
    [out_param userdat tmp res] = inputgui('title', guititle, 'geom', geometry, 'uilist',uilist,'eval', 'set(gcf,''tag'', ''pop_pacgui''); setappdata(findobj(''tag'', ''pop_pacgui''),''userdata'',struct(''freqindx1'', [], ''freqindx2'', []))', 'helpcom','pophelp(''pop_pac'');');
%%
    %%
    % Retreiving Inputs
    if isempty(res), return; end 
    
    indexfreqs1 = userdat.freqindx1;
    indexfreqs2 = userdat.freqindx2;
    pooldata    = datatypel_list{res.datatype};
    freqs1      = str2num(res.freq1);
    freqs2      = str2num(res.freq2);
    freqscale   = fastif(res.chckbx_logamp, 'log', 'linear');
    
    if numel(indexfreqs1) ~= numel(indexfreqs2)
        error('pop_pac: An equal number of components/channels mus be selected for getting the Phase and Amplitude data');
         return;
    end
      
    % Applying Bonferroni correction for the number of channel/components
    % computed
    if res.chckbx_stat
        if res.bonfcorr && ~isempty(res.pvalue)
            res.pvalue = str2num(res.pvalue) / (numel(indexfreqs1) * numel(indexfreqs2));
        else
            res.pvalue = str2num(res.pvalue);
        end
    else
        res.bonfcorr = 0;
        res.pvalue   = [];
    end
    
    % Updating g. options
    g.nboot     = str2num(res.nsurrogates);
    g.alpha     = res.pvalue;
    g.bonfcorr  = res.bonfcorr;
    g.method    = method_list{res.method};
    g.nfreqs1   = str2num(res.nfreqs1);
    g.nfreqs2   = str2num(res.nfreqs2)  ;
    g.freqscale = freqscale;
    if nsginstalled_flag
        g.compflag = fastif(res.chckbx_nsgt, 'nsg', 'local');
        tmpnsg = eval( [ '{' res.nsgopt '}' ] );
        if ~isempty( tmpnsg )
            for i = 1:2:numel(tmpnsg)
                if any(strcmp(tmpnsg{i},fieldnames(g)))
                    g.(tmpnsg{i}) = tmpnsg{i+1};
                end
            end
        end
    end
    
    % Options to call eeg_pac
    options    = {'freqs'     str2num(res.freq1)       'freqs2'     str2num(res.freq2)   'nfreqs1'   g.nfreqs1 ...
                  'nfreqs2'   g.nfreqs2                'freqscale'  g.freqscale          'method'     g.method...
                  'nboot'     g.nboot                  'alpha'      g.alpha              'bonfcorr'  g.bonfcorr'...
                  'tlimits'  minmax(EEG.times)};
    
    % Adding pop_pac options
    tmpparams = eval( [ '{' res.edit_optinput '}' ] );
    
    % Updating current function parameters and inputs to eeg_pac
    opttmp = tmpparams;
    if ~isempty( tmpparams )
        for i = 1:2:numel(opttmp)
            if any(strcmp(opttmp{i},fieldnames(g)))
                g.(opttmp{i}) = opttmp{i+1};
            else
                options{end+1} = opttmp{i};
                options{end+1} = opttmp{i+1};
            end
        end
    end
else
    options = {'freqs' freqs1 'freqs2' freqs2  'tlimits'  minmax(EEG.times)};
    options = {options{:} varargin{:}};
end


if strcmpi(g.compflag, 'local')
% Retreiving data
if strcmpi(pooldata,'channels')
    [dim_chan,~,dim_trial] = size(EEG.data); clear tmp;
    
    % Checking channel indices
    if ~all([all(ismember(indexfreqs1, 1:dim_chan)), all(ismember(indexfreqs2, 1:dim_chan))])
        error('Invalid data index');
    end
    
else % Component
    if isempty(EEG.icaact)
        EEG.icaact = eeg_getdatact( EEG,'component',1:size(EEG.icawinv,2));
    end
    [dim_comp,~,dim_trial] = size(EEG.data); clear tmp;
    
    % Checking component indices
    if ~all(all(ismember(indexfreqs1, 1:dim_comp)), all(ismember(indexfreqs2, 1:dim_comp)))
        error('Invalid data index');
    end
end

% Getting indices of the trials
if isempty(g.freq1_trialindx) || isempty(g.freq2_trialindx)
    g.freq1_trialindx = 1:dim_trial;
    g.freq2_trialindx = 1:dim_trial;
end

% Check if parameters are the same
%--
if ~g.cleanup
    methodindx  = find(strcmp(options,'method' ));
    opttmp = options;
    opttmp(methodindx:methodindx+1) = []; % removing method
    optionshash = gethashcode(std_serialize({'pooldata', pooldata, opttmp{:}}));
    if isfield(EEG.etc,'eegpac') && isfield(EEG.etc.eegpac, 'cache')
        [~, tmpstruct] = eeg_cache(EEG.etc.eegpac(1).cache, optionshash);
        if isempty(tmpstruct) 
            if ~g.forcecomp
                prompt = 'PACTools is about to remove previous saved PAC computations:';
                guititle = 'pop_pac - PAC results deletion';
                abortcomp = 'Exit';
                ignorewarn = 'Continue';
                answer = questdlg2(prompt,guititle,abortcomp,ignorewarn,ignorewarn);
                
                switch answer
                    case 'Continue'
                        disp('pop_pac: Proceeding to delete saved PAC results per user request');
                    case 'Exit'
                        disp('pop_pac: PAC computation has been canceled per user request. If you want to avoid this check use the ''forcecomp'' argument');
                        return;
                end
            end
            % Cleanup if different parameters
            EEG.etc.pacplotopt = [];
            EEG.etc.eegpac = [];
            EEG.etc.eegpac.cache = eeg_cache([], optionshash, {'pooldata', pooldata, opttmp{:}});
        end
    else
        EEG.etc.pacplotopt = [];
        EEG.etc.eegpac = [];
        EEG.etc.eegpac.cache = eeg_cache([], optionshash, {'pooldata', pooldata, opttmp{:}});
        
    end
else
    EEG.etc.eegpac = [];
    EEG.etc.pacplotopt = [];
end

if ~isempty(EEG.etc.eegpac) && isfield(EEG.etc.eegpac, 'dataindx')
      LastCell = length(EEG.etc.eegpac);
else, LastCell = 0;
end
c = LastCell+1;

for ichanpair = 1:length(indexfreqs1)
    
    % Check if measure is already computed
    ChanIndxExist = []; MeasureExist = 0;
    if ~isempty(EEG.etc.eegpac) && isfield(EEG.etc.eegpac, 'dataindx')
        ChanIndxExist = find(cell2mat(cellfun(@(x) all(x==[indexfreqs1(ichanpair) indexfreqs2(ichanpair)]), {EEG.etc.eegpac.dataindx}, 'UniformOutput', 0)));
        if ~isempty(ChanIndxExist)
            MeasureExist = isfield(EEG.etc.eegpac(ChanIndxExist),g.method) && ~isempty(EEG.etc.eegpac(ChanIndxExist).(g.method));
        end
    end
    
    if ~MeasureExist
        % Retreiving data for phase (X) and Amplitude
        if strcmpi(pooldata,'channels')
            X = squeeze(EEG.data(indexfreqs1(ichanpair),:,g.freq1_trialindx));
            Y = squeeze(EEG.data(indexfreqs2(ichanpair),:,g.freq2_trialindx));
            datatype = 1;
        else % Component
            X = squeeze(EEG.icaact(indexfreqs1(ichanpair),:,g.freq1_trialindx));
            Y = squeeze(EEG.icaact(indexfreqs2(ichanpair),:,g.freq2_trialindx));
            datatype = 2;
        end
        X = reshape(X,EEG.pnts,EEG.trials);
        Y = reshape(Y,EEG.pnts,EEG.trials);
        
        % Compute PAC
        [~, ~, freqs1, freqs2, alltfX, alltfY,~, pacstruct, tfXtimes, tfYtimes] = eeg_pac(X, Y, EEG.srate,  options{:});
        
        if ~isempty(ChanIndxExist)
                EEG.etc.eegpac(ChanIndxExist).(g.method) = pacstruct.(g.method);
        elseif ~isfield(EEG.etc.eegpac,'dataindx')
            EEG.etc.eegpac(1).(g.method) = pacstruct.(g.method);
            EEG.etc.eegpac(1).dataindx   = [indexfreqs1(ichanpair),indexfreqs2(ichanpair)] ;
            EEG.etc.eegpac(1).datatype   = datatype;
            EEG.etc.eegpac(1).params     = pacstruct.params;
            
            EEG.etc.eegpac(1).labels = {};
            if  EEG.etc.eegpac(1).datatype == 1
                if ~isempty(EEG.chanlocs)
                    EEG.etc.eegpac(1).labels = {[EEG.chanlocs(indexfreqs1(ichanpair)).labels '-' EEG.chanlocs(indexfreqs2(ichanpair)).labels]};
                else
                    EEG.etc.eegpac(1).labels = {['Chan' num2str(indexfreqs1(ichanpair)) '-' 'Chan' num2str(indexfreqs2(ichanpair))]};
                end
            else
                 EEG.etc.eegpac(1).labels = {['IC' num2str(indexfreqs1(ichanpair)) '-' 'IC' num2str(indexfreqs2(ichanpair))]};
            end
            c = c+1;
        else
            EEG.etc.eegpac(c).(g.method) = pacstruct.(g.method);
            EEG.etc.eegpac(c).dataindx   = [indexfreqs1(ichanpair),indexfreqs2(ichanpair)] ;
            EEG.etc.eegpac(c).datatype   = datatype;
            EEG.etc.eegpac(c).params     = pacstruct.params;
            EEG.etc.eegpac(c).cache      = EEG.etc.eegpac(1).cache;
            
            if  EEG.etc.eegpac(c).datatype == 1
                if ~isempty(EEG.chanlocs)
                    EEG.etc.eegpac(c).labels = {[EEG.chanlocs(indexfreqs1(ichanpair)).labels '-' EEG.chanlocs(indexfreqs2(ichanpair)).labels]};
                else
                    EEG.etc.eegpac(c).labels = {['Chan' num2str(indexfreqs1(ichanpair)) '-' 'Chan' num2str(indexfreqs2(ichanpair))]};
                end
            else
                EEG.etc.eegpac(c).labels = {['IC' num2str(indexfreqs1(ichanpair)) '-' 'IC' num2str(indexfreqs2(ichanpair))]};
            end
            c = c+1;
        end
    else
        disp(['pop_pac: Skipping computation of PAC for data pair with index: ' num2str(ichanpair) ]);
    end
end

% Creating fields for plottimg
EEG = pop_comodpacparams(EEG, 'default');
EEG = pop_comodtpacparams(EEG, 'default');
EEG = pop_tfpacparams(EEG, 'default');
EEG = pop_trialspacparams(EEG, 'default');

% Sorting fields and removing empty fields
% EEG = removeemptypac(EEG);
tmpval = setdiff( fieldnames(EEG.etc.eegpac),{'dataindx'    'datatype'  'labels'  'params' 'cache'})';
EEG.etc.eegpac = orderfields(EEG.etc.eegpac, {'dataindx'    'datatype'  'labels'  'params' 'cache' tmpval{:}});
com = sprintf('EEG = pop_pac(EEG,''%s'',[%s],[%s],[%s],[%s],%s);',pooldata,num2str([freqs1(1) freqs1(end)]),num2str([freqs2(1) freqs2(end)]),num2str(indexfreqs1),num2str(indexfreqs2), vararg2str(options(5:end)));
else
    %%%%%%%%%%%%%%%%%%%%%%
    %%% NSG Submission %%%
    %%%%%%%%%%%%%%%%%%%%%%
    try
        nsg_info;  % get information on where to create the temporary file
    catch
        error('Plugin nsgportal needs to be in the MATLAB path');
    end
    
    %  Section 1: Create temporary folder and save data
    tmpJobPath = fullfile(outputfolder, 'pactmp');
    if exist(tmpJobPath,'dir'), rmdir(tmpJobPath,'s'); end
    mkdir(tmpJobPath);
    
    % Save data in temporary folder previously created.
    % Here you may change the file name to match the one in the script you will run in NSG
    pop_saveset(EEG,'filename', EEG.filename, 'filepath', tmpJobPath);
    
    % Copy toolbox to folder. temporary until updated in NSG
    pactoolfolder = which('pop_pac.m');
    pacpath = fileparts(pactoolfolder);
    copyfile(pacpath,tmpJobPath);
    
    %  Section 2
    %  Manage m-file to be executed in NSG
    %  Write m-file to be run in NSG.
    %  Options defined in plugin are written into the file
    pac_singlesubj_writejobfile
    
    % Section 3
    % Submit job to NSG
    pop_nsg('run',tmpJobPath,'filename', 'pacssnsg_job.m', 'jobid', g.jobid,'runtime', g.runtime);
    display([char(10) 'PACTools job (jobID:'  g.jobid ') has been submitted to NSG' char(10) ...
                      'Copy or keep in mind the jobID assigned to this job to retreive the results later on.' char(10)...
                      'You may follow the status of your job through pop_nsg'...
        char(10)]);
    rmdir(tmpJobPath,'s');
    return;  
end
end

% function EEG = removeemptypac(EEG)
% fieldvals = fieldnames(EEG.etc.eegpac);
% for i =1: length(EEG.etc.eegpac)
%     emptyfield = find(structfun(@isempty, EEG.etc.eegpac{i}));
%     if ~isempty(emptyfield)
%         EEG.etc.eegpac{i} = rmfield(EEG.etc.eegpac({i}, fieldvals(emptyfield));
%     end 
% end
% end