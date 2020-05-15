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

STUDY = default_pacparams(STUDY);
TMPSTUDY = STUDY;
com = '';
if ~isfield(STUDY.etc, 'pac'), STUDY.etc.pac=[]; end
if isempty(varargin)
    
    rmsFlag = fastif(strcmpi(STUDY.etc.pac.tfpacparam.averagemode, 'rms'), 0, 1);
    icaFlag = fastif(isnan(STUDY.etc.pac.tfpacparam.topotime), 1, 0);
    
    tmpTimeRange = STUDY.etc.pac.tfpacparam.timerange;
    tmpFreqRange1 = STUDY.etc.pac.tfpacparam.freqrange1;
    tmpFreqRange2 = STUDY.etc.pac.tfpacparam.freqrange2;
    
    if strcmpi(STUDY.etc.pac.tfpacparam.averagechan,'off')
        multipleChansVal   = 1; % scalp array
    else
        multipleChansVal   = 2; % average channels
    end
    
    cb_multiplechan    = [ '    if ~isempty(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''timerange''), ''string'')))) ==1,' ...
                           '       set(findobj(gcbf, ''tag'', ''timerange''), ''string'', '''');' ...
                           '    end;' ...
                           '    if ~isempty(get(findobj(gcbf, ''tag'', ''freqrange1''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''freqrange1''), ''string'')))) ==1,' ...
                           '       set(findobj(gcbf, ''tag'', ''freqrange1''), ''string'', '''');' ...
                           '    end;' ...
                           '    if ~isempty(get(findobj(gcbf, ''tag'', ''freqrange2''), ''string'')) && length(unique(str2num(get(findobj(gcbf, ''tag'', ''freqrange2''), ''string'')))) ==1,' ...
                           '       set(findobj(gcbf, ''tag'', ''freqrange2''), ''string'', '''');' ...
                           '    end;'];    
    uilist = { ...
        {'style' 'text'       'string' 'Time-Freq PAC plotting options' 'fontweight' 'bold' 'tag', 'tfpac' } ...
        {'style' 'text'       'string' '(At least a single value must be defined for Phase or Amplitude frequencies)'  'tag', 'tfpac' } ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'}  {'style' 'edit'       'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Phase Freq. value or range in Hz [Low High]'} {'style' 'edit' 'string' num2str(tmpFreqRange1) 'tag' 'freqrange1' } ...
        {'style' 'text'       'string' 'Amp. Freq. value range in Hz [Low High]'} {'style' 'edit'  'string' num2str(tmpFreqRange2) 'tag' 'freqrange2' } ...
        {} ...
        {'style' 'text'       'string' 'Multiple channels selection' 'fontweight' 'bold' 'tag', 'spec' 'fontsize', 12} ...
        {} {'style' 'popupmenu'  'string' { 'Plot channels individually' 'Average of selected channels' } 'value' multipleChansVal 'tag' 'multiplechan' 'callback' cb_multiplechan } {} };
    evalstr = 'set(findobj(gcf, ''tag'', ''ersp''), ''fontsize'', 12);';
    otherline = [ 0.6 .4 ];
    chanline  = [ 0.07 0.8];
    geometry = { 1  1 otherline otherline otherline 1 1 chanline  1};
    
    if icaFlag
        uilist = uilist(1:end-2);
        geometry = geometry(1:end-2);
    end
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Time-Freq PAC plotting parameters -- pop_tfpacparams()', 'eval', evalstr);
    if isempty(res), return; end
    
    % decode input
    % ------------
    if ~isfield(res, 'multiplechan') res.multiplechan = 0; end
    res.timerange  = str2num( res.timerange );
    res.freqrange1 = str2num( res.freqrange1 );
    res.freqrange2 = str2num( res.freqrange2 );
    
    % build command call
    % ------------------
    options = {};

    if ~isequal(res.timerange, STUDY.etc.pac.tfpacparam.timerange),   options = { options{:} 'timerange' res.timerange };   end
    if ~isequal(res.freqrange1, STUDY.etc.pac.tfpacparam.freqrange1), options = { options{:} 'freqrange1' res.freqrange1 }; end
    if ~isequal(res.freqrange2, STUDY.etc.pac.tfpacparam.freqrange2), options = { options{:} 'freqrange1' res.freqrange2 }; end
    
    % mutliple channel option
    % -----------------------
    if res.multiplechan == 1
        if ~isequal('off', STUDY.etc.pac.tfpacparam.averagechan), options = { options{:} 'averagechan' 'off' }; end
    else
        if ~isequal('on', STUDY.etc.pac.tfpacparam.averagechan), options = { options{:} 'averagechan' 'on' }; end
    end
    
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_tfpacparams(STUDY, options{:});
        com = sprintf('STUDY = pop_tfpacparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_pacparams(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.pac.tfpacparam), 'exact'))
                STUDY.etc.pac.tfpacparam = setfield(STUDY.etc.pac.tfpacparam, varargin{index}, varargin{index+1});
            end
        end
    end
end

% scan clusters and channels to remove erspdata info if timerange etc. have
% changed (not neccesary here)

function STUDY = default_pacparams(STUDY)
    if ~isfield(STUDY.etc, 'pac'), STUDY.etc.pac = []; end
    if ~isfield(STUDY.etc.pac,'tfpacparam'), STUDY.etc.pac.tfpacparam = []; end
    if ~isfield(STUDY.etc.pac.tfpacparam, 'timerange'),     STUDY.etc.pac.tfpacparam.timerange = []; end
    if ~isfield(STUDY.etc.pac.tfpacparam, 'freqrange1'),    STUDY.etc.pac.tfpacparam.freqrange1 = []; end
    if ~isfield(STUDY.etc.pac.tfpacparam, 'freqrange2'),    STUDY.etc.pac.tfpacparam.freqrange2 = []; end
    if ~isfield(STUDY.etc.pac.tfpacparam, 'averagechan'),   STUDY.etc.pac.tfpacparam.averagechan = []; end