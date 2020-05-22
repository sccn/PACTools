% pop_comodtpacparams() - Set plotting and statistics parameters for 
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

function [ STUDY, com ] = pop_comodtpacparams(STUDY, varargin);

STUDY = default_comodtparams(STUDY);
com = '';
if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt=[]; end
if isempty(varargin)
      
    tmpTimeRange       = STUDY.etc.pacplotopt.comodtparam.timerange;
    tmpNPoints         = STUDY.etc.pacplotopt.comodtparam.npoints;
    tmpFreqRange1      = STUDY.etc.pacplotopt.comodtparam.freqrange1;
    tmpFreqRange2      = STUDY.etc.pacplotopt.comodtparam.freqrange2;
           
    uilist = { ...
        {'style' 'text'       'string' 'Comodulogram plotting options' 'fontweight' 'bold' 'tag', 'tfpac' } ...
        {'style' 'text'       'string' 'Time range ([Low High]) or value in ms'}  {'style' 'edit'  'string' num2str(tmpTimeRange) 'tag' 'timerange' } ...
        {'style' 'text'       'string' 'Number of latencies [Integer]'}  {'style' 'edit'  'string' num2str(tmpNPoints) 'tag' 'npoints' } ...
        {'style' 'text'       'string' 'Phase Freq. value or range in Hz [Low High]'} {'style' 'edit' 'string' num2str(tmpFreqRange1) 'tag' 'freqrange1' } ...
        {'style' 'text'       'string' 'Amp. Freq. value range in Hz [Low High]'} {'style' 'edit'  'string' num2str(tmpFreqRange2) 'tag' 'freqrange2' } {} };
    otherline = [ 0.6 .4 ];
    geometry = { 1 otherline otherline otherline otherline 1};
    
    [~, ~, ~, res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'skipline', 'off', 'title', 'Temporal Comodulogram plotting parameters -- pop_comodpacparams()');
    if isempty(res), return; end
    
    % decode input
    % ------------
    res.npoints       = str2num( res.npoints );
    res.timerange     = str2num( res.timerange );
    res.freqrange1    = str2num( res.freqrange1 );
    res.freqrange2    = str2num( res.freqrange2 );
    
    % build command call
    % ------------------
    options = {};
    if ~isequal(res.npoints,    STUDY.etc.pacplotopt.comodtparam.npoints),        options = { options{:} 'npoints'    res.npoints };    end
    if ~isequal(res.timerange,  STUDY.etc.pacplotopt.comodtparam.timerange),      options = { options{:} 'timerange'  res.timerange };  end
    if ~isequal(res.freqrange1, STUDY.etc.pacplotopt.comodtparam.freqrange1),     options = { options{:} 'freqrange1' res.freqrange1 }; end
    if ~isequal(res.freqrange2, STUDY.etc.pacplotopt.comodtparam.freqrange2),     options = { options{:} 'freqrange2' res.freqrange2 }; end
        
    % execute options
    % ---------------
    if ~isempty(options)
        STUDY = pop_comodtpacparams(STUDY, options{:});
        com = sprintf('STUDY = pop_comodtpacparams(STUDY, %s);', vararg2str( options ));
    end
else
    if strcmpi(varargin{1}, 'default')
        STUDY = default_comodtparams(STUDY);
    else
        for index = 1:2:length(varargin)
            if ~isempty(strmatch(varargin{index}, fieldnames(STUDY.etc.pacplotopt.comodtparam), 'exact'))
                STUDY.etc.pacplotopt.comodtparam = setfield(STUDY.etc.pacplotopt.comodtparam, varargin{index}, varargin{index+1});
            end
        end
    end
end

% scan clusters and channels to remove erspdata info if timerange etc. have
% changed (not neccesary here)

function STUDY = default_comodtparams(STUDY)
    if ~isfield(STUDY.etc, 'pacplotopt'), STUDY.etc.pacplotopt = []; end
    if ~isfield(STUDY.etc.pacplotopt,'comodtparam'), STUDY.etc.pacplotopt.comodtparam = []; end
    if ~isfield(STUDY.etc.pacplotopt.comodtparam, 'npoints'),           STUDY.etc.pacplotopt.comodtparam.npoints      = 10; end
    if ~isfield(STUDY.etc.pacplotopt.comodtparam, 'timerange'),         STUDY.etc.pacplotopt.comodtparam.timerange    = []; end
    if ~isfield(STUDY.etc.pacplotopt.comodtparam, 'freqrange1'),        STUDY.etc.pacplotopt.comodtparam.freqrange1   = []; end
    if ~isfield(STUDY.etc.pacplotopt.comodtparam, 'freqrange2'),        STUDY.etc.pacplotopt.comodtparam.freqrange2   = []; end