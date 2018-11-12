% pop_pac() - Compute phase-amplitude coupling.
%
% Usage:
%   >>  pac = pop_pac(EEG);
%   >>  pac = pop_pac(EEG,)
%
% Inputs:
%  EEG         - input dataset
%  pooldata    - ('channels' || 'component' ) 
%  phasefreq   -
%  ampfreq     -
%  indexphase  -
%  indexamp    -
%
% Optional inputs:
%
%
%
% Outputs:
%
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2016
%
% Copyright (C) 2016  Ramon Martinez-Cancino,INC, SCCN
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

function [EEG,com] = pop_pac(EEG,pooldata,phasefreq,ampfreq,indexphase,indexamp,varargin)
com = [];
flag_compute = true;

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
    disp('std_infocluster() error: calling convention {''key'', value, ... } error'); return;
end
try g.freq1_trialindx;      catch, g.freq1_trialindx  =  1:dim_trial; end % Opt only for pop_pac
try g.freq2_trialindx;      catch, g.freq2_trialindx  =  1:dim_trial; end % Opt only for pop_pac

try g.nboot;                catch, g.nboot            =  200;         end
try g.alpha;                catch, g.alpha            =  0.05;        end
try g.bonfcorr;             catch, g.bonfcorr         =  0;           end
try g.methodpac;            catch, g.methodpac        =  'glm';       end
try g.nfreqs1;              catch, g.nfreqs1          =  1;           end
try g.nfreqs2;              catch, g.nfreqs2          =  1;           end
try g.cleanup,              catch, g.cleanup          = 0;            end
try g.freqscale,            catch, g.freqscale        = 'linear';     end % may be 'linear' or 'log'
try g.ptspercent,          catch, g.ptspercent       = 0.05;         end


if nargin < 6 
    % Defaults for GUI
    phasefreq      = [4 15];
    ampfreq        = [EEG.srate/2-20 EEG.srate/2-2];
    nfreqs1        = length(phasefreq(1):2:phasefreq(2));
    nfreqs2        = length(ampfreq(1):3:ampfreq(2)) ;
    freq1_dataindx = 1;
    freq2_dataindx = 1;
    
    % Scan here for add-in methods and add to gui list
    
    method_defaultlistgui = {'Mean vector length modulation index (Canolty et al.)',...
                             'Kullback-Leibler modulation index (Tort et al.)',...
                             'General linear model (Penny et al.)',...
                             'Instantaneous MIPAC', 'Event related MIPAC'};
    method_defaultlist    = {'mvlmi','klmi','glm','instmipac', 'ermipac'};
    % ---
    % Checking if ICA
    if isempty(EEG.icawinv)
        datatypel_list = {'Channels'} ;
    else
        datatypel_list = {'Components','Channels'} ;
    end
    
    % Here define other types of CFC
    % Note: In case of adding more methods, add the modality here at the
    % end of the cell array, and then update 'data1_list' and 'data2_list'
    % in the callback 'callback_setdata'
    
    cfctype_list = {'Phase-Amp.', 'Amp.-Amp.', 'Phase-Phase'}; 
 
   
    method_listgui = {'Mean vector length modulation index (Canolty et al.)', 'Kullback-Leibler modulation index (Tort et al.)','General linear model (Penny et al.)', 'Instantaneous MIPAC', 'Event related MIPAC'};
    method_list    = {'mvlmi','klmi','glm','instmipac', 'ermipac'};
    guititle = 'pop_pac() - Test for event-related Phase/Amplitude Coupling (PAC)'; 
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

    callback_setdata = ['data1_list = {''Phase Data'' ''Amp. Data'' ''Phase Data''};' ... 
                        'data2_list = {''Amp. Amp.''  ''Amp. Data'' ''Phase Data''};' ...
                         'tmpval = get(findobj(''tag'', ''datatyppac''),''value'');'...
                         'set(findobj(''tag'', ''data1''), ''string'',data1_list{tmpval});'...
                         'set(findobj(''tag'', ''data2''), ''string'',data2_list{tmpval})'];
    
    geometry = { {2.7 9 [0 0]    [1 1]}  {2.7 9 [0.3  0] [0.8 1] }  {2.7 9 [1.4 0]  [1 1]}    {2.7 9 [1.75 0] [0.8 1] }...
                 {2.7 9 [0.6 1]  [1 1]}  {2.7 9 [1.25 1] [1 1]}     {2.7 9 [1.9 1]  [1 1]} ...
                 {2.7 9 [0.18 2] [1 1]}  {2.7 9 [0.6  2] [0.65 1]}  {2.7 9 [1.25 2] [0.65 1]} {2.7 9 [1.9 2] [0.6 1]}...
                 {2.7 9 [0.18 3] [1 1]}  {2.7 9 [0.6  3] [0.65 1]}  {2.7 9 [1.25 3] [0.65 1]} {2.7 9 [1.9 3] [0.6 1]}...
                 {2.7 9 [0 4]    [1 1]}  {2.7 9 [0.58 4] [1.96 1]}...
                 {2.7 9 [0 5]    [1 1]}  {2.7 9 [1.25 5] [1.25 1]} ...
                 {2.7 9 [0 6]    [1 1]}  {2.7 9 [0.58 6] [1.96 1]} {2.7 9 [0.7  6] [1.96 1]}...
                 {2.7 9 [0.2 7]  [1 1]}  {2.7 9 [1.24 7] [0.5 1]}  {2.7 9 [1.7  7] [1 1]} {2.7 9 [2  7] [0.5 1]}...
                 {2.7 9 [0.2 8]  [1 1]}  {2.7 9 [1.24 8] [0.5 1]}...
                 {2.7 9 [0.2 9]  [1 1]}  {2.7 9 [1.24 9] [0.5 1]}...
                 {2.7 9 [0  10]  [1 1]}};
 %{2.7 9 [1.26 7] [0.5 1]} {2.7 9 [1.26 7] [0.5 1]}  {'style' 'edit' 'string' 1/num2str(g.ptspercent) 'tag' 'nblog_edit' 'enable' 'off'}
    
    uilist = {{'style' 'text' 'string' 'Data type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' datatypel_list 'tag' 'datatype' 'value' 1} {'style' 'text' 'string' 'CFC type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' cfctype_list 'tag' 'datatyppac' 'value' 1 'callback' callback_setdata}...% 
              {'style' 'text' 'string' 'Comp/chan indices' 'fontweight' 'normal'} {'style' 'text' 'string' 'Freq range [lo hi] Hz' 'fontweight' 'normal'}         {'style' 'text' 'string' '   # Frequencies' 'fontweight' 'normal'} ...
              {'style' 'text' 'string' 'Phase data' 'fontweight' 'normal' 'tag' 'data1'} ...
                                                                                                {'style' 'edit' 'string' num2str(freq1_dataindx) 'tag' 'freq1_dataindx'}...
                                                                                                {'style' 'edit' 'string' num2str(phasefreq)      'tag' 'freq1'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs1)        'tag' 'nfreqs1'}...
              {'style' 'text' 'string' 'Amp data ' 'fontweight' 'normal' 'tag' 'data2'} ...
                                                                                                {'style' 'edit' 'string' num2str(freq2_dataindx) 'tag' 'freq2_dataindx'}...
                                                                                                {'style' 'edit' 'string' num2str(ampfreq)        'tag' 'freq2'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs2)        'tag' 'nfreqs2'}...
              {'style' 'text' 'string' 'PAC method' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' method_listgui 'tag' 'method'}...
               {'style' 'text' 'string' 'Comand line options' 'fontweight' 'bold'} {'style' 'edit' 'string' ' ' 'tag' 'edit_optinput'}...
              {'style' 'text' 'string' 'PAC statistics' 'fontweight' 'bold'} {'style' 'checkbox' 'tag' 'chckbx_stat' 'callback' callback_chkcbxstat 'value' 0} {'style' 'text' 'string' '(off)' 'fontweight' 'normal' 'tag' 'label_statstate'}...
              {'style' 'text' 'string' '# Surrogates' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.nboot) 'tag' 'nsurrogates' 'enable' 'off'}  {'style' 'text' 'string' '# blocks' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(1/g.ptspercent) 'tag' 'nblocks_edit' 'enable' 'off'} ...
              {'style' 'text' 'string' 'Significance threshold (0<p<1)' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.alpha) 'tag' 'pvalue' 'enable' 'off'}... 
              {'style' 'text' 'string' 'Correct for multiple comparisons' 'callback' 'close(gcbf);' } {'style' 'checkbox' 'tag' 'bonfcorr' 'enable' 'off'}...
              {}};
          
    [out_param userdat tmp res] = inputgui('title', guititle, 'geom', geometry, 'uilist',uilist, 'helpcom','pophelp(''pop_pac'');');
%%
    % Retreiving Inputs
    if isempty(res), return; end 
    
    indexphase = str2num(res.freq1_dataindx);
    indexamp   = str2num(res.freq2_dataindx);
    pooldata   = datatypel_list{res.datatype};
    phasefreq  = str2num(res.freq1);
    ampfreq    = str2num(res.freq2);
    
    % Applying Bonferroni correction for the number of channel/components
    % computed
    if res.chckbx_stat
        if res.bonfcorr && ~isempty(res.pvalue)
            res.pvalue = str2num(res.pvalue) / (numel(indexphase) * numel(indexamp));
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
    g.methodpac =  method_list{res.method};
    g.nfreqs1   = str2num(res.nfreqs1);
    g.nfreqs2   = str2num(res.nfreqs2)  ;
    
    % Options to call eeg_pac
    options    = {'freqs' str2num(res.freq1)          'freqs2'    str2num(res.freq2)          'methodpac' method_list{res.method}...
                  'nboot' str2num(res.nsurrogates)    'alpha'     res.pvalue,                 'nfreqs1',str2num(res.nfreqs1), ...
                  'nfreqs2',str2num(res.nfreqs2)      'bonfcorr'  res.bonfcorr'};
    
    % Adding options
    tmpparams = eval( [ '{' res.edit_optinput '}' ] );
    options   = {options{:} tmpparams{:}};
else
    options = {'freqs' phasefreq 'freqs2' ampfreq};
    options = {options{:} varargin{:}};
end

% Retreiving data
if strcmpi(pooldata,'channels')
    [dim_chan,~,dim_trial] = size(EEG.data); clear tmp;
    
    % Checking channel indices
    if ~all(all(ismember(indexphase, 1:dim_chan)), all(ismember(indexamp, 1:dim_chan)))
        error('Invalid data index');
    end
    
    else % Component
    if isempty(EEG.icaact)
        EEG.icaact = eeg_getdatact( EEG,'component',1:length(EEG.icawinv));
    end
    [dim_comp,~,dim_trial] = size(EEG.data); clear tmp;
    
    % Checking component indices
    if ~all(all(ismember(indexphase, 1:dim_comp)), all(ismember(indexamp, 1:dim_comp)))
        error('Invalid data index');
    end
end

% Getting indices of the trials
if isempty(g.freq1_trialindx) || isempty(g.freq2_trialindx)
    g.freq1_trialindx = 1:dim_trial;
    g.freq2_trialindx = 1:dim_trial;
end

% Check if this was already computed (not implemented)

%--

timefreq_phase = cell(1,length(indexphase));
timefreq_amp   = cell(1,length(indexamp));

if ~isfield(EEG.etc,'eegpac') || g.cleanup
    EEG.etc.eegpac = [];
elseif isfield(EEG.etc.eegpac,'params')
    
    if strcmpi(g.freqscale, 'log')      
        inputparams.freqs_phase = linspace(log(phasefreq(1)), log(phasefreq(end)), g.nfreqs1);
        inputparams.freqs_phase = exp(inputparams.freqs_phase);
        
        inputparams.freqs_amp = linspace(log(ampfreq(1)), log(ampfreq(end)), g.nfreqs2);
        inputparams.freqs_amp = exp(inputparams.freqs_amp);        
    else
        inputparams.freqs_phase = linspace(phasefreq(1), phasefreq(2), g.nfreqs1); % this should be OK for FFT
        inputparams.freqs_amp   = linspace(ampfreq(1), ampfreq(2), g.nfreqs2);
    end
    
     inputparams.srate             = EEG.srate;
     inputparams.signif.alpha      = g.alpha;
     inputparams.signif.nboot      = g.nboot;
     inputparams.signif.bonfcorr   = g.bonfcorr;
     
     tmpstrc = EEG.etc.eegpac.params;
     tmpstrc.signif = rmfield(tmpstrc.signif,'ptspercent');
    if ~isequal(EEG.etc.eegpac.params, inputparams)
        EEG.etc.eegpac = [];
    end  
end

c = 1;
for ichan_phase = 1:length(indexphase)
    % Retreiving data for phase (X)
    if strcmpi(pooldata,'channels')
        X = EEG.data(indexphase(ichan_phase),:,g.freq1_trialindx);
    else % Component
        X = EEG.icaact(indexphase(ichan_phase),:,g.freq1_trialindx);
    end
    %----    
    for ichan_amp = 1:length(indexamp)
        % Retreiving data for amplitude (Y)
        if strcmpi(pooldata,'channels')
            Y = EEG.data(indexamp(ichan_amp),:,g.freq2_trialindx);
        else % Component
            Y = EEG.icaact(indexamp(ichan_amp),:,g.freq2_trialindx);
        end
        %----       
        % Running eeg_pac
        % Three options, so we can save the TF decompositions us them later in the loop
        if all([isempty(timefreq_phase{ichan_phase}),isempty(timefreq_amp{ichan_amp})])
            [~, timesout, freqs1, freqs2, alltfX, alltfY,~, pacstruct] = eeg_pac(X', Y', EEG.srate,  options{:});
            
            % populating  phase cell
            timefreq_phase{ichan_phase}.timesout = timesout;
            timefreq_phase{ichan_phase}.freqs    = freqs1;
            timefreq_phase{ichan_phase}.alltf    = alltfX;
            
            timefreq_amp{ichan_amp}.timesout = timesout;
            timefreq_amp{ichan_amp}.freqs    = freqs2;
            timefreq_amp{ichan_amp}.alltf    = alltfY;
            
        elseif isempty(timefreq_amp{ichan_amp})
            tmpopt = [options 'alltfXstr' timefreq_phase{ichan_phase}] ;
            [~, timesout, ~, freqs2, ~, alltfY,~, pacstruct] = eeg_pac(X', Y', EEG.srate,  tmpopt{:}); clear tmpopt;
            
            % populating  amp cell
            timefreq_amp{ichan_amp}.timesout = timesout;
            timefreq_amp{ichan_amp}.freqs    = freqs2;
            timefreq_amp{ichan_amp}.alltf    = alltfY;
            
        elseif isempty(timefreq_phase{ichan_phase})
            tmpopt = [options 'alltfYstr' timefreq_amp{ichan_amp}] ;  
            [~, timesout, freqs1, ~,alltfX, ~,~, pacstruct] = eeg_pac(X', Y', EEG.srate,  tmpopt{:}); clear tmpopt;
            % populating  phase cell
            timefreq_phase{ichan_phase}.timesout = timesout;
            timefreq_phase{ichan_phase}.freqs    = freqs1;
            timefreq_phase{ichan_phase}.alltf    = alltfX;
        end
         EEG.etc.eegpac.(g.methodpac){c} = pacstruct.(g.methodpac);       
         EEG.etc.eegpac.freqcell{c}   = [ichan_phase,ichan_amp] ; 
         c = c+1;
    end
end

% Common stuff
EEG.etc.eegpac.params = pacstruct.params;
     
com = sprintf('pop_pac(EEG,''%s'',[%s],[%s],[%s],[%s],%s);',pooldata,num2str(phasefreq),num2str(ampfreq),num2str(indexphase),num2str(indexamp), vararg2str(options(5:end)));
end