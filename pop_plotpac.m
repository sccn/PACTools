% pop_plotpac() - Call GUI to compute cross-frequency-coupling coupling.
%             Second level function to compute CFC by calling eeg_pac.m
% Usage:
%   >>  pac = pop_plotpac(EEG);
%
% Inputs:
%  EEG          - [Structure] Input dataset as an EEGLAB EEG structure

% Outputs:
%  h            - Handles of figures generated
%
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

function [EEG,com] = pop_plotpac(EEG, plottype, pacmethod, phasedataindx, ampdataindx, varargin)


try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end;
catch
    disp('pop_plotpac() error: calling convention {''key'', value, ... } error'); return;
end;

try g.plotopt;                                 catch, g.plotopt         = {};              end;
try g.plotsignif;                              catch, g.plotsignif      = 0;               end;
try g.fampval;                                 catch, g.fampval         = [];              end;
try g.fphaseval;                               catch, g.fphaseval       = [];              end;
try g.timeval;                                 catch, g.timeval         = [];              end;

% Check if EEG file 
if ~isfield(EEG,'etc') && ~isfield(EEG,'data')
    error('pop_plotpac() error: Ivalid EEG structure provided as input');
end
% Check if PAC was computed
if ~isfield(EEG.etc,'eegpac') && isempty(EEG.etc.eegpac)
    error('pop_plotpac() error: PAC must be computed');
end

if nargin < 6 
 %%   
fieldnamesval = fieldnames(EEG.etc.eegpac);
if length(fieldnamesval)>3
    listmethod   = fieldnamesval(4:end); % getting computed method
else
    error('pop_plotpac() error: PAC must be computed');
end
listchanindx     = cellfun(@(x) cellstr(num2str(x)),EEG.etc.eegpac.dataindx)';


listplot   = {'Comodulogram',...
              'PhaseAmpTime','Amp-PhaseTime','Phase-AmpTime','Time-PhaseAmp',...
              'AmpPhase-TrialTime','Amp-PhaseTrialTime','Phase-AmpTrialTime' };
          
% Building the GUI
    guititle = 'pop_plotpac()';
    uilist = { ...
        { 'style', 'text', 'string','Method','fontweight' 'bold'},{'Style', 'popupmenu', 'string', listmethod , 'tag', 'listbox_method'},{ 'style', 'text', 'string','Chan/Comp. Indx','fontweight' 'bold'},{'Style', 'popupmenu', 'string', listchanindx , 'tag', 'listbox_chanindx'},...
        { 'style', 'text', 'string','Plot type','fontweight' 'bold'}, {'Style', 'popupmenu', 'string', listplot , 'tag', 'listbox_plot'}...
        { 'style', 'text', 'string','Phase freq. (Hz)','fontweight' 'bold'},{ 'Style', 'edit', 'string','','tag', 'edit_phase'}...
        { 'style', 'text', 'string','Amp freq. (Hz)','fontweight' 'bold'},  { 'Style', 'edit', 'string','','tag', 'edit_amp'},...
        { 'style', 'text', 'string','Time (s)','fontweight' 'bold'},      { 'Style', 'edit', 'string','','tag', 'edit_time'},...
        { 'style', 'text', 'string','Command line options','fontweight' 'bold'},{ 'Style', 'edit', 'string','','tag', 'edit_opt'},...
         };
 
    geometry = {{4 4 [0 0]    [1 1]}  {4 4 [0.6 0] [1.2 1] }  {4 4 [2 0]  [0.8 1]}  {4 4 [2.8 0] [1.2 1] }...
                {4 4 [0 1]    [1 1]}  {4 4 [0.6 1] [1.2 1]}...
                {4 4 [0 2]    [1 1]}  {4 4 [0.9 2] [0.5 1]}...
                {4 4 [1.5 2]  [1 1]}  {4 4 [2.3 2] [0.5 1]}...
                {4 4 [2.9 2]  [1 1]}  {4 4 [3.4 2] [0.5 1]}...
                {4 4 [0 3]    [2 1]}  {4 4 [1.4 3] [2.5 1]}};

            [out_param userdat tmp res] = inputgui('title', guititle, 'geom', geometry, 'uilist',uilist, 'helpcom','pophelp(''pop_plotpac'');');
            
            plottype      = listplot{res.listbox_plot};
            pacmethod     = listmethod{res.listbox_method};
            tmpchanindx   = str2num(listchanindx{res.listbox_chanindx});
            phasedataindx = tmpchanindx(1);
            ampdataindx   = tmpchanindx(2);
            g.fampval       = str2num(res.edit_amp);
            g.fphaseval     = str2num(res.edit_phase);
            g.timeval       = str2num(res.edit_time);
            
            if ~isempty(res.edit_opt)
                tmpparams = eval( [ '{' res.edit_opt '}' ] );
                g.plotopt   = {tmpparams{:}};
            end
end
  
fighandle = eeg_plotpac(EEG,plottype,'pacmethod',pacmethod,...
                                     'phasedataindx',phasedataindx,...
                                     'ampdataindx',ampdataindx,...
                                     'plotopt', g.plotopt,...
                                     'plotsignif', g.plotsignif,...
                                     'fampval', g.fampval,...
                                     'fphaseval', g.fphaseval,...
                                     'timeval', g.timeval);