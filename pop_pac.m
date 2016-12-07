% pop_pac() - Compute phase-amplitude coupling.
%
% Usage:
%   >>  pac = pop_pac(EEG);
%   >>  pac = pop_pac(EEG,)
%
% Inputs:
%
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
if isempty(EEG)
    fprintf(2,'pop_pac error : Empty input provided \n');
    return;
end
EEG = eeg_checkset(EEG);

% get inputs here
[tmp,tmp,dim_trial] = size(EEG.data); clear tmp;

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('std_infocluster() error: calling convention {''key'', value, ... } error'); return;
end;
try g.freq1_trialindx;      catch, g.freq1_trialindx  =  1:dim_trial; end;
try g.freq2_trialindx;      catch, g.freq2_trialindx  =  1:dim_trial; end;
try g.nboot;                catch, g.nboot            =  200;           end;
try g.alpha;                catch, g.alpha            =  0.05;        end;
try g.methodpac;            catch, g.methodpac        =  'mvlmi';     end;
try g.forcecomp;            catch, g.forcecomp        =  1;           end;
try g.nfreqs1;              catch, g.nfreqs1          =  1;           end;
try g.nfreqs2;              catch, g.nfreqs2          =  1;           end;

if nargin <8 
    % Defaults for GUI
    phasefreq      = [4 15];
    ampfreq        = [20 EEG.srate/2];
    nfreqs1        = length(phasefreq(1):2:phasefreq(2));
    nfreqs2        = length(ampfreq(1):3:ampfreq(2)) ;
    freq1_dataindx = [1 2];
    freq2_dataindx = [1 2];
    
    % Checking if ICA
    if isempty(EEG.icawinv)
        datatypel_list = {'Channels'} ;
    else
        datatypel_list = {'Components','Channels'} ;
    end
    method_listgui = {'Mean vector length modulation index (Canolty et al.)', 'Kullback-Leibler modulation index (Tort et al.)','Sin/cos regression (Penny et al.)'};
    method_list    = {'mvlmi','klmi','glm'};
    guititle = 'pop_pac() - Test for event-related Phase/Amplitude Coupling (PAC)'; 
    geometry = {    {2.7 8 [0 0] [1 1]}  {2.7 8 [0.6 0] [0.8 1]}...
                                    {2.7 8 [0.6 1] [1 1]}    {2.7 8 [1.25 1] [1 1]}    {2.7 8 [1.9 1] [1 1]} ...
                {2.7 8 [0.18 2] [1 1]} {2.7 8 [0.6 2] [0.65 1]} {2.7 8 [1.25 2] [0.65 1]} {2.7 8 [1.9 2] [0.6 1]}...
                {2.7 8 [0.18 3] [1 1]} {2.7 8 [0.6 3] [0.65 1]} {2.7 8 [1.25 3] [0.65 1]} {2.7 8 [1.9 3] [0.6 1]}...
                {2.7 8 [0 4] [1 1]} {2.7 8 [0.58 4] [1.96 1]}...
                {2.7 8 [0 5] [1 1]} ...
                {2.7 8 [0.2 6] [1 1]} {2.7 8 [1.26 6] [0.5 1]}...
                {2.7 8 [0.2 7] [1 1]} {2.7 8 [1.26 7] [0.5 1]}};
 %%
    
    uilist = {{'style' 'text' 'string' 'PAC data type' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' datatypel_list 'tag' 'datatype' 'value' 1}...
              {'style' 'text' 'string' 'Comp/chan indices' 'fontweight' 'normal'} {'style' 'text' 'string' 'Freq range [lo hi] Hz' 'fontweight' 'normal'}         {'style' 'text' 'string' '   # Frequencies' 'fontweight' 'normal'} ...
              {'style' 'text' 'string' 'Phase data' 'fontweight' 'normal'} ...
                                                                                                {'style' 'edit' 'string' num2str(freq1_dataindx) 'tag' 'freq1_dataindx'}...
                                                                                                {'style' 'edit' 'string' num2str(phasefreq)      'tag' 'freq1'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs1)        'tag' 'nfreqs1'}...
              {'style' 'text' 'string' 'Amp data ' 'fontweight' 'normal'} ...
                                                                                                {'style' 'edit' 'string' num2str(freq2_dataindx) 'tag' 'freq2_dataindx'}...
                                                                                                {'style' 'edit' 'string' num2str(ampfreq)        'tag' 'freq2'}...
                                                                                                {'style' 'edit' 'string' num2str(nfreqs2)        'tag' 'nfreqs2'}...
              {'style' 'text' 'string' 'PAC method' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' method_listgui 'tag' 'method'}...
              {'style' 'text' 'string' 'PAC statistics' 'fontweight' 'bold'} ...
              {'style' 'text' 'string' '# Surrogates' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.nboot) 'tag' 'nsurrogates'}...
              {'style' 'text' 'string' 'Significance threshold (0<p<1)' 'callback' 'close(gcbf);' } {'style' 'edit' 'string' num2str(g.alpha) 'tag' 'pvalue'}};
    
    [out_param userdat tmp res] = inputgui('title', guititle, 'geom', geometry, 'uilist',uilist, 'helpcom','pophelp(''pop_pac'');');
   
    % Retreiving Inputs
    if isempty(res), return; end;
    options    = {'freqs' str2num(res.freq1)          'freqs2' str2num(res.freq2)          'methodpac' method_list{res.method}...
                  'nboot' str2num(res.nsurrogates)    'alpha' str2num(res.pvalue),         'nfreqs1',str2num(res.nfreqs1), ...
                  'nfreqs2',str2num(res.nfreqs2)      };
              
    indexphase = str2num(res.freq1_dataindx);
    indexamp   = str2num(res.freq2_dataindx);
    pooldata   = datatypel_list{res.datatype};
    phasefreq  = str2num(res.freq1);
    ampfreq    = str2num(res.freq2);
else
    options = {'freqs' phasefreq 'freqs2' ampfreq 'nfreqs1' g.nfreqs1 'nfreqs2' g.nfreqs2 'methodpac' g.methodpac 'nboot' g.nboot 'alpha' g.alpha };
end

% Retreiving data
if strcmpi(pooldata,'channels')
    [dim_chan,tmp,dim_trial] = size(EEG.data); clear tmp;
    
    % Checking channel indices
    if ~all(all(ismember(indexphase, 1:dim_chan)), all(ismember(indexamp, 1:dim_chan)))
        error('Invalid data index');
    end
    
    else % Component
    if isempty(EEG.icaact)
        EEG.icaact = eeg_getdatact( EEG,'component',1:length(EEG.icawinv));
    end
    [dim_comp,tmp,dim_trial] = size(EEG.data); clear tmp;
    
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

% Check if this was already computed
optextended = {pooldata, phasefreq, ampfreq, indexphase, indexamp, options,varargin};

if  ~g.forcecomp && isfield(EEG.etc,'eegpac') && ~isempty(EEG.etc.eegpac)
    serialized_inputs = std_serialize(optextended);
    computedflag = isequal(serialized_inputs,EEG.etc.eegpac.serialize);
    if computedflag
        flag_compute = false;
        EEG.etc.eegpac.serialize = serialized_inputs;
    end;
end

if flag_compute
     timefreq_phase = cell(1,length(indexphase));
     timefreq_amp   = cell(1,length(indexamp));
     
    for ichan_phase = 1:length(indexphase)
        % Retreiving data for phase (X)
        if strcmpi(pooldata,'channels')         
            X = EEG.data(indexphase(ichan_phase),:,g.freq1_trialindx);
        else % Component             
            X = EEG.icaact(indexphase(ichan_phase),:,g.freq1_trialindx);
        end
        
        for ichan_amp = 1:length(indexamp)
            % Retreiving data for amplitude (Y)
            if strcmpi(pooldata,'channels')
                Y = EEG.data(indexamp(ichan_amp),:,g.freq2_trialindx);
            else % Component
                Y = EEG.icaact(indexamp(ichan_amp),:,g.freq2_trialindx);
            end
            
            % Running eeg_pac
            % Four options, so we can save the TF decompositions us them later in the loop
            if all([isempty(timefreq_phase{ichan_phase}),isempty(timefreq_amp{ichan_amp})])
                [pacval, timesout, freqs1, freqs2, tmp, alltfX, alltfY,crossfcoh_pval, pacstruct] = eeg_pac(X, Y, EEG.srate,  options{:}); clear tmp;
                
                % populating  phase cell
                timefreq_phase{ichan_phase}.timesout = timesout;
                timefreq_phase{ichan_phase}.freqs    = freqs1;
                timefreq_phase{ichan_phase}.alltf    = alltfX;
                
                timefreq_amp{ichan_amp}.timesout = timesout;
                timefreq_amp{ichan_amp}.freqs    = freqs2;
                timefreq_amp{ichan_amp}.alltf    = alltfY;
                
            elseif isempty(timefreq_amp{ichan_amp})
                tmpopt = options;
                tmpopt{end+1} = 'alltfXstr' ;
                tmpopt{end+1} = timefreq_phase{ichan_phase};
                
                [pacval, timesout, freqs1, freqs2, tmp, alltfX, alltfY,crossfcoh_pval, pacstruct] = eeg_pac(X, Y, EEG.srate,  tmpopt{:}); clear tmp;
                
                % populating  amp cell
                timefreq_amp{ichan_amp}.timesout = timesout;
                timefreq_amp{ichan_amp}.freqs    = freqs2;
                timefreq_amp{ichan_amp}.alltf    = alltfY;
                
            elseif isempty(timefreq_phase{ichan_phase})
                tmpopt = options;
                tmpopt{end+1} = 'alltfYstr' ;
                tmpopt{end+1} = timefreq_amp{ichan_amp};
             
                [pacval, timesout, freqs1, freqs2, tmp, alltfX, alltfY,crossfcoh_pval, pacstruct] = eeg_pac(X, Y, EEG.srate,  tmpopt{:}); clear tmp;
                % populating  phase cell
                timefreq_phase{ichan_phase}.timesout = timesout;
                timefreq_phase{ichan_phase}.freqs    = freqs1;
                timefreq_phase{ichan_phase}.alltf    = alltfX;
            else
                tmpopt = options;
                tmpopt{end+1} = 'alltfXstr' ;
                tmpopt{end+1} = timefreq_phase{ichan_phase};
                tmpopt{end+1} = 'alltfYstr' ;
                tmpopt{end+1} = timefreq_amp{ichan_amp};
             
                [pacval, timesout, freqs1, freqs2, tmp, alltfX, alltfY,crossfcoh_pval, pacstruct] = eeg_pac(X, Y, EEG.srate,  tmpopt{:}); clear tmp;
            end
            
            % Storing results in EEG
            if isfield(EEG.etc,'eegpac') && all([ichan_phase==1,ichan_amp ==1])
                EEG.etc.eegpac = [];
            end
            
%             pacstruct.indxchanphase = indexphase(ichan_phase);
%             pacstruct.indxchanamp   = indexamp(ichan_amp);
%             pacstruct.options       = options;
%             pacstruct.pacval        = pacval;
%             pacstruct.timesout      = timesout;
%             pacstruct.freqs1        = freqs1;
%             pacstruct.freqs2        = freqs2;
            
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.pacstruct      = indexphase; % Deprecate after finish development
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.indxchanphase  = indexphase(ichan_phase);
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.indxchanamp    = indexamp(ichan_amp);
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.pacval         = pacval;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.pval           = pacstruct.pval;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.signifmask     = pacstruct.signifmask;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.peakangle      = pacstruct.peakangle;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.betas          = pacstruct.betas;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.bin_average    = pacstruct.bin_average;
            EEG.etc.eegpac.datapac{ichan_phase,ichan_amp}.composites     = pacstruct.composites;
        end
    end
    
    % Common stuff
    EEG.etc.eegpac.pahseindx    = indexphase;
    EEG.etc.eegpac.ampindx      = indexamp;
    
    EEG.etc.eegpac.options      = options;
    EEG.etc.eegpac.method       = pacstruct.method;
    EEG.etc.eegpac.freqs_phase  = pacstruct.freqs_phase;
    EEG.etc.eegpac.freqs_amp    = pacstruct.freqs_amp;
    EEG.etc.eegpac.alpha        = pacstruct.alpha;
    EEG.etc.eegpac.ptspercent   = pacstruct.ptspercent;
    EEG.etc.eegpac.nboots       = pacstruct.nboots;
    EEG.etc.eegpac.srate        = pacstruct.srate;
    EEG.etc.eegpac.timesout     = timesout;
    EEG.etc.eegpac.options      = options;
    EEG.etc.eegpac.freqs1       = freqs1;
    EEG.etc.eegpac.freqs2       = freqs2;
    EEG.etc.eegpac.nbinskl      = pacstruct.nbinskl;
          
end

com = sprintf('pop_pac(EEG,''%s'',[%s],[%s],[%s],[%s],%s);',pooldata,num2str(phasefreq),num2str(ampfreq),num2str(indexphase),num2str(indexamp), vararg2str(options(5:end)));
end