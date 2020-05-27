% pop_tfpacparams() - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean (and optionally 
%                    single-trial) ERSP and ITC measures and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erspparams) which is used
%                    whenever plotting is performed by the function
%                    std_erspplot.
% Usage:    
%   >> STUDY = pop_tfpacparams(STUDY, 'key', 'val', ...);   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% ERSP/ITC image plotting options:
%  'timerange'   - [min max] ERSP/ITC plotting latency range in ms. 
%                  {default: the whole output latency range}.
%  'freqrange1'  - [min max] PAC plotting frequency range in ms. 
%                  {default: the whole output frequency range}
%  'freqrange2'  - [min max] ERSP/ITC plotting frequency range in ms. 
%                  {default: the whole output frequency range}
%  'averagechan' - ['on'|'off'] average/rms data channels when several are
%                  selected ('on') or plot them individually ('off')
%
% See also:  std_erspplot()
%
% Authors: Ramon Martinez-Cancino
%          Arnaud Delorme

% Copyright (C) Ramon Martinez-Cancino, 2020
%
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ STUDY, com ] = pop_tfpacparams(STUDY, varargin);

STUDY = default_tfpacparams(STUDY);
com = '';
if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt=[]; end
if isempty(varargin)
      
    tmpTFixFreq        = STUDY.etc.pacplotopt.tfparam.fixfreq;
    tmpTimeRange       = STUDY.etc.pacplotopt.tfparam.timerange;
    tmpFreqRange       = STUDY.etc.pacplotopt.tfparam.freqrange;
    tmpFreqVal1        = STUDY.etc.pacplotopt.tfparam.freqval1;
    tmpFreqVal2        = STUDY.etc.pacplotopt.tfparam.freqval2;
    
    if tmpTFixFreq==1, chbx_freq1 = 1, chbx_freq2 = 0; enable_chbx1 =  'on'; enable_chbx2 = 'off';
    else,              chbx_freq1 = 0, chbx_freq2 = 1; enable_chbx1 = 'off'; enable_chbx2 = 'on';
    end
    cb_chbox1 = 'set(findobj(''tag'',''freqval2''),''enable'' ,''off''); set(findobj(''tag'',''freqval1''),''enable'' ,''on'');set(findobj(''tag'',''fix_freq2''),''value'' ,0);'; 
    cb_chbox2 = 'set(findobj(''tag'',''freqval1''),''enable'' ,''off''); set(findobj(''tag'',''freqval2''),''enable'' ,''on'');set(findobj(''tag'',''fix_freq1''),''value'' ,0);';
    
    uilist = { ...
        {'style' 'text'       'string' 'PAC Freq-Time plotting options' 'fontweight' 'bold' 'tag', 'tfpac' } ...
        {'style' 'text'       'string' 'Define frequency value' 'fontweight' 'bold'}  {'style' 'text' 'string' 'Fix value' 'fontweight' 'bold'}...
        {'style' 'text'       'string' 'Phase Freq.'}  {'style' 'edit' 'string' num2str(tmpFreqVal1) 'tag' 'freqval1' 'enable' enable_chbx1}   {'style' 'checkbox' 'value' chbx_freq1 'tag' 'fix_freq1'  'callback' cb_chbox1} ...
        {'style' 'text'       'string' 'Amplitude Freq.'}  {'style' 'edit' 'string' num2str(tmpFreqVal2) 'tag' 'freqval2'  'enable' enable_chbx2}   {'style' 'checkbox' 'value' chbx_freq2 'tag' 'fix_freq2' 'callback' cb_chbox2}...
        {'style' 'text'       'string' 'Time range in ms [Low High]'}  {'style' 'edit'       'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Frequency range in Hz [Low High]'} {'style' 'edit' 'string' num2str(tmpFreqRange) 'tag' 'freqrange' } {} };     
    otherline = [ 0.2 .2 ];
    chbxline = [ .3 .2 .1 ];
    geometry = { 1  [.7 .2] chbxline chbxline otherline otherline 1 };
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Time-frequency PAC plotting parameters -- pop_tfpacparams()');
                                        
    if isempty(res), return; end
    if isempty(res.freqval1) && isempty(res.freqval2),  disp('A value must be defined for Amplitude or Phase frequency'); return; end
    
    % decode input
    % ------------
    res.timerange          = str2num(res.timerange) ;
    res.freqrange          = str2num(res.freqrange) ;
    if res.fix_freq1 == 1,  res.fixfreq = 1; res.freqval1 = str2num(res.freqval1) ; else, res.freqval1 = []; end
    if res.fix_freq2 == 1,  res.fixfreq = 2; res.freqval2 = str2num(res.freqval2) ; else, res.freqval2 = []; end
    
    % build command call
    % ------------------
    options = {};
    
    if ~isequal(res.timerange, STUDY.etc.pacplotopt.tfparam.timerange),       options = { options{:} 'timerange'  res.timerange };  end
    if ~isequal(res.freqrange, STUDY.etc.pacplotopt.tfparam.freqrange),       options = { options{:} 'freqrange'  res.freqrange };  end
    if ~isequal(res.freqval1,  STUDY.etc.pacplotopt.tfparam.freqval1),        options = { options{:} 'freqval1'   res.freqval1 };   end
    if ~isequal(res.freqval2,  STUDY.etc.pacplotopt.tfparam.freqval2),        options = { options{:} 'freqval2'   res.freqval2 };   end
    if ~isequal(res.fixfreq,   STUDY.etc.pacplotopt.tfparam.fixfreq),         options = { options{:} 'fixfreq'    res.fixfreq };    end
        
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_tfpacparams(STUDY, options{:});
        if isstudy(STUDY)
            structname = 'STUDY';
        else
            structname = 'EEG';
        end
        com = sprintf('%s = pop_tfpacparams(%s, %s);',structname, structname, vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_tfpacparams(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.pacplotopt.tfparam), 'exact'))
                STUDY.etc.pacplotopt.tfparam = setfield(STUDY.etc.pacplotopt.tfparam, varargin{index}, varargin{index+1});
            end
        end
    end
end

function STUDY = default_tfpacparams(STUDY)

phval = []; % Fill out when implementing STUDY level
if ~isstudy(STUDY)
    try
    phval = ceil(numel(STUDY.etc.eegpac(1).params.freqs_phase)/2);
    catch
        disp('Input must be valid a EEG/STUDY structure');
        return;
    end
end    
if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt = []; end
if ~isfield(STUDY.etc.pacplotopt,'tfparam'), STUDY.etc.pacplotopt.tfparam = []; end
if ~isfield(STUDY.etc.pacplotopt.tfparam, 'timerange'),          STUDY.etc.pacplotopt.tfparam.timerange  = [];    end
if ~isfield(STUDY.etc.pacplotopt.tfparam, 'freqrange'),          STUDY.etc.pacplotopt.tfparam.freqrange  = [];    end
if ~isfield(STUDY.etc.pacplotopt.tfparam, 'freqval1'),           STUDY.etc.pacplotopt.tfparam.freqval1   = phval; end
if ~isfield(STUDY.etc.pacplotopt.tfparam, 'freqval2'),           STUDY.etc.pacplotopt.tfparam.freqval2   = [];    end
if ~isfield(STUDY.etc.pacplotopt.tfparam, 'fixfreq'),            STUDY.etc.pacplotopt.tfparam.fixfreq    = 1;     end