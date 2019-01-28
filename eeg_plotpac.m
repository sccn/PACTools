% eeg_plotpac() - Plots for phase amplitude coupling data generated
%                       with ERPAC Tool
%
% Usage:
%   >> eeg_plotpac(EEG,plottype);
%   >> eeg_plotpac(EEG,plottype,'key1', 'val1', 'key2', val2' ...);
%
% Inputs:
%    EEG       - (structure) EEG structure or pac structure(see specification on toolbox README.md)
%    plottype  - (string)    Plot type. This argument should be one of the following plot types 
%
%                 For 1D data: 
%                 'Comodulogram': Standard comodulogram Phase Vs Amplitude
% 
%                 For 3D data(single trial with time):
%                 'PhaseAmpTime': Tri-dimensional plot of Phase Vs Amplitude Vs Time
%                 'Amp-PhaseTime': Bi-dimensional plot of Phase Vs Time for a pinned Amplitude value
%                 'Phase-AmpTime': Bi-dimensional plot of Amplitude Vs Time for a pinned Phase value
%                 'Time-PhaseAmp': Bi-dimensional plot of Phase Vs Amp (comodulogram) for a pinned Time value.
% 
%                 For 4D data (Multiple trials with time)
%                 'AmpPhase-TrialTime': 
%                 'Amp-PhaseTrialTime':
%                 'Phase-AmpTrialTime':
%    
%                 For customized plots:
%                               
%   Optional inputs
%     plotopt
%     plotsignif
%     cell2plot
%     phasechanindx
%     ampchanindx
%     pacmethod
%     fampval
%     fphaseval
%     timeval

% Outputs:
%         h   - Handles of the figures generated
%                       
% Author: Ramon Martinez Cancino
%
% See also:
%
% Copyright (C) 2019 Ramon Martinez Cancino, UCSD, INC, SCCN
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


function hfig = eeg_plotpac(EEG,plottype, varargin)
hfig = [];

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end;
catch
    disp('eeg_plotpac() error: calling convention {''key'', value, ... } error'); return;
end;

try g.plotopt;                                 catch, g.plotopt         = {};              end;
try g.plotsignif;                              catch, g.plotsignif      = 0;               end;
try g.cell2plot;                               catch, g.cell2plot       = 1;               end;
try g.phasechanindx;                           catch, g.phasechanindx   = [];              end;
try g.ampchanindx;                             catch, g.ampchanindx     = [];              end;
try g.pacmethod;                               catch, g.pacmethod       = [];              end;

try g.fampval;                                 catch, g.fampval         = [];              end;
try g.fphaseval;                               catch, g.fphaseval       = [];              end;
try g.timeval;                                 catch, g.timeval         = [];              end;

% Detect if EEG or pacstruct structure
if isfield(EEG,'etc') && isfield(EEG,'data')
    if isfield(EEG.etc,'eegpac') && ~isempty(EEG.etc.eegpac)  
        % Reconstruct pacstruct from EEG structure
        pacstruct.params = EEG.etc.eegpac.params;       
        % if not pacmethod provided use the first available
        if isempty(g.pacmethod)
            fieldnameftmp = fieldnames(EEG.etc.eegpac);
            g.pacmethod = fieldnameftmp{4};
        end
        methodstructtmp = EEG.etc.eegpac.(g.pacmethod);
        % Determining  cell to plot
        if ~isempty(g.phasechanindx) && ~isempty(g.ampchanindx)
            g.cell2plot = find(cell2mat(cellfun(@(x) isequal(x,[g.phasechanindx,g.ampchanindx]), EEG.etc.eegpac.chanindx, 'UniformOutput', 0)));
        end
        pacstruct.(g.pacmethod) = methodstructtmp{g.cell2plot};  
    else
        error('eeg_plotpac() error: Missing or empty structure ''eegpac''');
    end
else
    pacstruct = EEG;
end

% List of plots
plotlistdim1 = {'Comodulogram'};
plotlistdim2 = {'PhaseAmpTime', 'Amp-PhaseTime', 'Phase-AmpTime', 'Time-PhaseAmp'};
plotlistdim3 = {'AmpPhase-TrialTime', 'Amp-PhaseTrialTime', 'Phase-AmpTrialTime'};
      
% Checkinf if plottype provided is in the list of valid plots
if ~ismember(plottype, {plotlistdim1{:} plotlistdim2{:} plotlistdim3{:}})
     error('eeg_plotpac() error: Invalid value for input ''plottype''');
end

% Checking if type of plot requested fits the dimensions 
invalidplot_errorflag = 0;
if (pacstruct.(g.pacmethod).dim == 1 && ~ismember(plottype,plotlistdim1))  
 invalidplot_errorflag = 1;
elseif pacstruct.(g.pacmethod).dim == 2 && ~ismember(plottype,plotlistdim2)
    invalidplot_errorflag = 1;
elseif pacstruct.(g.pacmethod).dim == 3 && ~ismember(plottype,plotlistdim3)
    invalidplot_errorflag = 1;
end
if invalidplot_errorflag
    error(['eeg_plotpac() error: Invalid ''plottype'' for PAC dimension equal to ' num2str(pacstruct.(g.pacmethod).dim)]);
end

 % Plotting start here
 switch  plottype
     % Comodulogram
     case 'Comodulogram'
         if g.plotsignif
             if ~isempty(pacstruct.(g.pacmethod).signif.signifmask)
                 g.plotopt(end+1:end+2) = {'signifmask', pacstruct.(g.pacmethod).signif.signifmask};     
             else
                 disp('eegplot pac() message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
             end
         end 
         
         if length(pacstruct.params.freqs_phase)==1 || length(pacstruct.params.freqs_amp)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         plotopttmp = g.plotopt;
         g.plotopt = {'title', 'Modulation Index'}; 
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         hfig = comodulogram(pacstruct.params.freqs_phase,pacstruct.params.freqs_amp,pacstruct.(g.pacmethod).pacval, g.plotopt{:});
     
     case 'PhaseAmpTime'
         if length(pacstruct.params.freqs_phase)==1 || length(pacstruct.params.freqs_amp)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         plotopttmp = g.plotopt;
         g.plotopt(end+1:end+10) = {'title', 'Modulation Index', 'zlabel', 'Phase Frequency (Hz)','ylabel', 'Amplitude Frequency (Hz)', 'xlabel', 'Time(sec)', 'cbartext', 'Modulation Index'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         hfig = comodulogramt(pacstruct.params.freqs_phase,pacstruct.params.freqs_amp,pacstruct.(g.pacmethod).times, pacstruct.(g.pacmethod).pacval, g.plotopt{:});
      
     case 'Amp-PhaseTime'   
         if length(pacstruct.params.freqs_phase)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         if isempty(g.fampval)
             fampvalindx = round(length(pacstruct.params.freqs_amp)/2);
         else
             % Check if input is valid
             if g.fampval<pacstruct.params.freqs_amp(1) || g.fampval>pacstruct.params.freqs_amp(end)
                 error('eeg_plotpac(): Invalid value for input ''fampval''');
             end
             [trash,fampvalindx]= min(abs(pacstruct.params.freqs_amp-g.fampval));
         end
         data = squeeze(pacstruct.(g.pacmethod).pacval(:,fampvalindx,:));
         times = pacstruct.(g.pacmethod).times;
         titl  = ['Phase Vs Time (f_{Amp} = ' num2str(pacstruct.params.freqs_amp(fampvalindx)) 'Hz)'];
         plotopttmp = g.plotopt;
         g.plotopt = {'erp', 'on'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         
         figure; 
         [trash,trash,trash,trash,axhndls] = erpimage(data',[],times,titl,0,0,'img_trialax_label','Amplitude Frequency (Hz)',...
                                              'cbar', 'on',...
                                              'cbar_title', 'Mod. Index',...
                                              'erp', 'on',...
                                              'yerplabel', 'Marginal MI', g.plotopt{:}) ;        
         
         set(axhndls{1},'YTickLabel',pacstruct.params.freqs_phase(get(axhndls{1},'YTick')));
         set(axhndls{1}, 'box', 'on'); set(axhndls{2}, 'box', 'on'); set(axhndls{3}, 'box', 'on');
         
     case 'Phase-AmpTime'
         if length(pacstruct.params.freqs_amp)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         if isempty(g.fphaseval)
             fphasevalindx = round(length(pacstruct.params.freqs_phase)/2);
         else
             % Check if input is valid
             if g.fphaseval<pacstruct.params.freqs_phase(1) || g.fphaseval>pacstruct.params.freqs_phase(end)
                 error('eeg_plotpac(): Invalid value for input ''fphaseval''');
             end
             [trash,fphasevalindx]= min(abs(pacstruct.params.freqs_phase-g.fphaseval));
         end
         data = squeeze(pacstruct.(g.pacmethod).pacval(fphasevalindx,:,:));
         times = pacstruct.(g.pacmethod).times;
         titl  = ['Phase Vs Time (f_{hase} = ' num2str(pacstruct.params.freqs_phase(fphasevalindx)) 'Hz)'];
         plotopttmp = g.plotopt;
         g.plotopt = {'erp', 'on'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         figure; 
         [trash,trash,trash,trash,axhndls] = erpimage(data',[],times,titl,0,0,'img_trialax_label','Phase Frequency (Hz)',...
                                              'cbar', 'on',...
                                              'cbar_title', 'Mod. Index',...
                                              'erp', 'on',...
                                              'yerplabel', 'Marginal MI', g.plotopt{:}) ;        
         
         set(axhndls{1},'YTickLabel',pacstruct.params.freqs_amp(get(axhndls{1},'YTick')));
         set(axhndls{1}, 'box', 'on'); set(axhndls{2}, 'box', 'on'); set(axhndls{3}, 'box', 'on');
         
     case 'Time-PhaseAmp'
         if length(pacstruct.params.freqs_phase)==1 || length(pacstruct.params.freqs_amp)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         if isempty(g.timeval)
             timeindx = 1;
         else
             % Check if input is valid
             minmaxval = minmax(pacstruct.(g.pacmethod).times);
             if g.timeval<minmaxval(1) || g.timeval>minmaxval(2)
                 error(['eeg_plotpac(): Invalid value for input ''timeval''. Value out of valid range (' num2str(minmaxval(1)) ',' num2str(minmaxval(2)) ')']);
             end
             % Finding timeindx
             [trash, timeindx] = min(abs(pacstruct.(g.pacmethod).times-g.timeval));
         end
         
         if g.plotsignif
             if ~isempty(pacstruct.(g.pacmethod).signif.signifmask)
                 g.plotopt(end+1:end+2) = {'signifmask', pacstruct.(g.pacmethod).signif.signifmask(:,:,timeindx)};     
             else
                 disp('eegplot pac() message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
             end
         end 
         
         plotopttmp = g.plotopt;
         g.plotopt = {'title', ['Modulation Index (time = ' num2str(pacstruct.(g.pacmethod).times(timeindx)) 's)']};    
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         hfig = comodulogram(pacstruct.params.freqs_phase,pacstruct.params.freqs_amp,squeeze(pacstruct.(g.pacmethod).pacval(:,:,timeindx)), g.plotopt{:});
         
     % AmpPhase-TrialTime     
     case 'AmpPhase-TrialTime'
         % Amp
         if isempty(g.fampval)
             fampvalindx = round(length(pacstruct.params.freqs_amp)/2);
         else
             % Check if input is valid
             if g.fampval<pacstruct.params.freqs_amp(1) || g.fampval>pacstruct.params.freqs_amp(end)
                 error('eeg_plotpac(): Invalid value for input ''fampval''');
             end
             [trash,fampvalindx]= min(abs(pacstruct.params.freqs_amp-g.fampval));
         end
         
         % Phase
         if isempty(g.fphaseval)
             fphasevalindx = round(length(pacstruct.params.freqs_phase)/2);
         else
             % Check if input is valid
             if g.fphaseval<pacstruct.params.freqs_phase(1) || g.fphaseval>pacstruct.params.freqs_phase(end)
                 error('eeg_plotpac(): Invalid value for input ''fphaseval''');
             end
             [trash,fphasevalindx]= min(abs(pacstruct.params.freqs_phase-g.fphaseval));
         end
         
         data = squeeze(pacstruct.(g.pacmethod).pacval(fphasevalindx,fampvalindx,:,:));
         times = pacstruct.(g.pacmethod).times;
         titl  = ['Trials Vs Time (f_{Phase} = ' num2str(pacstruct.params.freqs_phase(fphasevalindx)) 'Hz, f_{Amp} = ' num2str(pacstruct.params.freqs_amp(fampvalindx)) 'Hz )'];
         
         plotopttmp = g.plotopt;
         g.plotopt = {'erp', 'on'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         figure; 
         [trash,trash,trash,trash,axhndls] = erpimage(data',[],times,titl,0,0,'img_trialax_label','Trials',...
                                              'cbar', 'on',...
                                              'cbar_title', 'Mod. Index',...
                                              'erp', 'on',...
                                              'yerplabel', 'Marginal MI', g.plotopt{:}) ;        
         
         set(axhndls{1}, 'box', 'on'); set(axhndls{2}, 'box', 'on'); set(axhndls{3}, 'box', 'on');
     
     % Amp-PhaseTrialTime    
     case 'Amp-PhaseTrialTime'
         if length(pacstruct.params.freqs_phase)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         
         % Amp
         if isempty(g.fampval)
             fampvalindx = round(length(pacstruct.params.freqs_amp)/2);
         else
             % Check if input is valid
             if ~ismember(g.fampval,pacstruct.params.freqs_amp)
                 error('eeg_plotpac(): Invalid value for input ''fampval''');
             end
         end
         
         plotopttmp = g.plotopt;
         titl = ['Modulation Index (f_{Amp} ='  num2str(pacstruct.params.freqs_amp(fampvalindx))  'Hz)'];
         g.plotopt = {'title', titl, 'zlabel', 'Phase Frequency (Hz)','ylabel', 'Trials', 'xlabel', 'Time(sec)', 'cbartext', 'Modulation Index'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         hfig = comodulogramt(pacstruct.params.freqs_phase,pacstruct.params.freqs_amp,pacstruct.(g.pacmethod).times, squeeze(pacstruct.(g.pacmethod).pacval(:,fampvalindx,:,:)), g.plotopt{:});
     
     % Phase-AmpTrialTime
     case 'Phase-AmpTrialTime'
          if length(pacstruct.params.freqs_amp)==1
             error('eeg_plotpac() error: Insuficient number of frequencies for plotting');
         end
         % Phase
         if isempty(g.fphaseval)
             fphasevalindx = round(length(pacstruct.params.freqs_phase)/2);
         else
             % Check if input is valid
             if ~ismember(g.fphaseval,pacstruct.params.freqs_phase)
                 error('eeg_plotpac(): Invalid value for input ''fphaseval''');
             end
         end
         
         plotopttmp = g.plotopt;
         titl = ['Modulation Index (f_{Phase} ='  num2str(pacstruct.params.freqs_phase(fphasevalindx))  'Hz)'];
         g.plotopt = {'title', titl, 'zlabel', 'Amp. Frequency (Hz)','ylabel', 'Trials', 'xlabel', 'Time(sec)', 'cbartext', 'Modulation Index'};
         g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
         hfig = comodulogramt(pacstruct.params.freqs_phase,pacstruct.params.freqs_amp,pacstruct.(g.pacmethod).times, squeeze(pacstruct.(g.pacmethod).pacval(fphasevalindx,:,:,:)), g.plotopt{:});
 end
end