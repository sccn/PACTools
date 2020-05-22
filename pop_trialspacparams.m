% pop_trialspacparams - Set plotting and statistics parameters for 
%                    computing and plotting STUDY mean (and optionally 
%                    single-trial) ERSP and ITC measures and measure 
%                    statistics. Settings are stored within the STUDY 
%                    structure (STUDY.etc.erspparams) which is used
%                    whenever plotting is performed by the function
%                    std_erspplot.
% Usage:    
%   >> STUDY = pop_comodpacparams(STUDY, 'key', 'val', ...);   
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

function [ STUDY, com ] = pop_trialspacparams(STUDY, varargin)

STUDY = default_trialspacparams(STUDY);
com = '';
if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt=[]; end
if isempty(varargin)
      
    tmpTrialRange    = STUDY.etc.pacplotopt.trialparam.trialrange;
    tmpTimeRange     = STUDY.etc.pacplotopt.trialparam.timerange;
    tmpFreqVal1      = STUDY.etc.pacplotopt.trialparam.freqval1;
    tmpFreqVal2      = STUDY.etc.pacplotopt.trialparam.freqval2;
    
    uilist = { ...
        {'style' 'text'       'string' 'Trial-based PAC plotting options' 'fontweight' 'bold'} ...
        {'style' 'text'       'string' 'Phase Freq. value in Hz'}  {'style' 'edit' 'string' num2str(tmpFreqVal1) 'tag' 'freqval1' } ...
        {'style' 'text'       'string' 'Amp. Freq. value in Hz'}   {'style' 'edit' 'string' num2str(tmpFreqVal2) 'tag' 'freqval2' } ...
        {'style' 'text'       'string' 'Time range in ms [Low High]'}   {'style' 'edit' 'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Trial range'}  {'style' 'edit' 'string' tmpTrialRange 'tag' 'trialrange' } {} };
    otherline = [ 0.6 .4 ];
    geometry = { 1 otherline otherline otherline otherline 1};
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', ...
                                            'title', 'Comodulogram plotting parameters -- pop_comodpacparams()');
    if isempty(res), return; end
    
    % decode input
    % ------------
    res.timerange           = str2num( res.timerange );
    res.trialrange          = str2num( res.trialrange );
    res.freqval1            = str2num( res.freqval1 );
    res.freqval2            = str2num( res.freqval2 );
    
    % build command call
    % ------------------
    options = {};
    if ~isequal(res.trialrange,     STUDY.etc.pacplotopt.trialparam.trialrange),       options = { options{:} 'trialrange'      res.trialrange };    end
    if ~isequal(res.timerange,      STUDY.etc.pacplotopt.trialparam.timerange),        options = { options{:} 'timerange'       res.timerange };     end
    if ~isequal(res.freqval1,       STUDY.etc.pacplotopt.trialparam.freqval1),         options = { options{:} 'freqval1'        res.freqval1 };      end
    if ~isequal(res.freqval2,       STUDY.etc.pacplotopt.trialparam.freqval2),         options = { options{:} 'freqval2'        res.freqval2 };      end
        
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_trialspacparams(STUDY, options{:});
        if isstudy(STUDY)
            structname = 'STUDY';
        else
            structname = 'EEG';
        end
        com = sprintf('%s = pop_trialspacparams(%s, %s);',structname, structname, vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_trialspacparams(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.pacplotopt.trialparam), 'exact'))
                STUDY.etc.pacplotopt.trialparam = setfield(STUDY.etc.pacplotopt.trialparam, varargin{index}, varargin{index+1});
            end
        end
    end
end

function STUDY = default_trialspacparams(STUDY)
    if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt = []; end
    if ~isfield(STUDY.etc.pacplotopt,'trialparam'), STUDY.etc.pacplotopt.trialparam = []; end
  
    if ~isfield(STUDY.etc.pacplotopt.trialparam, 'trialrange'),         STUDY.etc.pacplotopt.trialparam.trialrange   = []; end
    if ~isfield(STUDY.etc.pacplotopt.trialparam, 'timerange'),          STUDY.etc.pacplotopt.trialparam.timerange    = []; end
    if ~isfield(STUDY.etc.pacplotopt.trialparam, 'freqval1'),           STUDY.etc.pacplotopt.trialparam.freqval1     = []; end
    if ~isfield(STUDY.etc.pacplotopt.trialparam, 'freqval2'),           STUDY.etc.pacplotopt.trialparam.freqval2     = []; end