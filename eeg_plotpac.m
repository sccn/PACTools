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
%                 'PhaseAmpTime':  Tri-dimensional plot of Phase Vs Amplitude Vs Time
%                 'Amp-PhaseTime': Bi-dimensional plot of Phase Vs Time for a pinned Amplitude value
%                 'Phase-AmpTime': Bi-dimensional plot of Amplitude Vs Time for a pinned Phase value
%                 'Time-PhaseAmp': Bi-dimensional plot of Phase Vs Amp (comodulogram) for a pinned Time value.
% 
%                 For 4D data (Multiple trials with time)
%                 'AmpPhase-TrialTime': Trials Vs Time for pinned values of Aplitude and Phase
%                 'Amp-PhaseTrialTime': Tridimensional plot Phase Vs Trials Vs Time for pinned value of Amplitude
%                 'Phase-AmpTrialTime': Tridimensional plot Amplitude Vs Trials Vs Times for pinned Phase values
%    
%                 For customized plots:
%                               
%   Optional inputs
%     plotopt
%     plotsignif
%     cell2plot
%     phasedataindx
%     ampdataindx
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
try g.phasedataindx;                           catch, g.phasedataindx   = [];              end;
try g.ampdataindx;                             catch, g.ampdataindx     = [];              end;
try g.pacmethod;                               catch, g.pacmethod       = [];              end;

try g.fampval;                                 catch, g.fampval         = [];              end;
try g.fphaseval;                               catch, g.fphaseval       = [];              end;
try g.timeval;                                 catch, g.timeval         = [];              end;

% Detect if EEG or pacstruct structure
if isfield(EEG,'etc') && isfield(EEG,'data')
    if isfield(EEG.etc,'eegpac') && ~isempty(EEG.etc.eegpac)  
        % Reconstruct pacstruct from EEG structure
        pacstruct.params = EEG.etc.eegpac(1).params; 
        
        % Determining  cell to plot
        if ~isempty(g.phasedataindx) && ~isempty(g.ampdataindx)
            g.cell2plot = find(cell2mat(cellfun(@(x) isequal(x,[g.phasedataindx,g.ampdataindx]), {EEG.etc.eegpac.dataindx}, 'UniformOutput', 0)));
        end
        
        % if not pacmethod provided use the first available
        if isempty(g.pacmethod)
            fieldnameftmp = fieldnames(EEG.etc.eegpac);
            tmpval = find(structfun(@(x) ~isempty(x),EEG.etc.eegpac(g.cell2plot)));
            tmpval2 = find(tmpval>4);
            g.pacmethod = fieldnameftmp{tmpval(tmpval2)};
        end
        pacstruct.(g.pacmethod) = EEG.etc.eegpac(g.cell2plot).(g.pacmethod);
          
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
         if ~isempty(g.timeval)
             g.plotopt(end+1:end+2) = {'timeval', g.timeval};
         end
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

% Auxiliar Functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function h = comodulogram(freqs_phase,freqs_amp,pacval,varargin)
    icadefs;
    h = [];
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('eeg_visualize_pac() error: calling convention {''key'', value, ... } error'); return;
end;

try g.signifmask; catch, g.signifmask          = [];              end;
try g.alphadata;  catch, g.alphadata           = 0.5;              end;
try g.title;      catch, g.title               = 'Modulation Index';              end;
    
    h = figure('Name', 'Modulation Index','Units','Normalized','Position', [3.4740e-01 2.4352e-01 3.4219e-01 5.4074e-01], 'Resize', 'off','Tag','comod_plot');
%     h_axes1 = axes('parent',h, 'Units', 'Normalized',...
%                    'ActivePositionProperty', 'outerposition',...
%                    'Color','None',...
%                    'Visible', 'off');
%     haxes = axes('parent',h,'Units', 'Normalized',...
%                 'Position',h_axes1.Position,...
%                 'ActivePositionProperty', 'outerposition');
    

    haxes = axes('parent',h,'Units', 'Normalized',...
                'ActivePositionProperty', 'outerposition');
            
    imagesc(freqs_phase,freqs_amp,pacval, 'parent', haxes);
    set(haxes,'YDir','normal');
    title(haxes, sprintf(g.title));
    colorbar(haxes);
    xlabel(haxes,'Phase Frequency [Hz]');
    ylabel(haxes,'Amplitude Frequency [Hz]');
    set(get(h,'Children'),'Fontsize',AXES_FONTSIZE_L+5);

     % Plot significance mask 
    if ~isempty(g.signifmask)
        haxes2 = axes('Position',haxes.Position,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
        haxes2.ActivePositionProperty = 'outerposition';
        linkaxes([haxes,haxes2]);
        h_trans = imagesc(freqs_phase,freqs_amp,g.signifmask', 'parent', haxes2);
        set(h_trans,'AlphaData',g.alphadata);
        set(haxes2,'Color','None','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','')
%         colormap(haxes,'jet');
        colormap(haxes2,'gray');
    end
    box on; grid on;
end
   
function h = comodulogramt(freqs_phase,freqs_amp,timevect, pacval, varargin)

icadefs;
h = [];
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('comodulogramt() error: calling convention {''key'', value, ... } error'); return;
end;

try g.ntimepoints;                     catch, g.ntimepoints = 10;                     end;
try g.timeval;                         catch, g.timeval = [];                           end;
try g.title;                           catch, g.title = 'Modulation Index';           end;
try g.zlabel;                          catch, g.zlabel = 'Phase Frequency (Hz)';      end;
try g.ylabel;                          catch, g.ylabel = 'Amplitude Frequency (Hz)';  end;
try g.xlabel;                          catch, g.xlabel = 'Time(sec)';                 end;
try g.cbartext;                        catch, g.cbartext = 'Modulation Index';        end; 
try g.comodtazimuth;                   catch, g.comodtazimuth = -10;                  end;
try g.comodtelevation;                 catch, g.comodtelevation = 23;                 end;

% Determining nearest time
if isempty(g.timeval)
    times = min(timevect):((max(timevect)-min(timevect))/(g.ntimepoints-1)):max(timevect);
else
    times = g.timeval;
end
[~, timeidx] = min(abs(repmat(timevect,length(times),1)-repmat(times',1,length(timevect))),'',2);
tplot = timevect(timeidx);   

%% Plotting start here

h    = figure('Name', g.title,'Units','Normalized','Position', [0.2349    0.3093    0.5547    0.2907],'Tag','comodt_plot');
haxes = axes('Units', 'Normalized','Color','None','parent', h);

% Plot
Z = pacval(:,:,timeidx);
[M,N,P] = size(Z);
for i=1:P
    % Create a plane at x=i
    hsurf = surface(tplot(i)*ones(1,M),1:N,repmat([M:-1:1],N,1),repmat([M:-1:1],N,1));
    % set the color of the plane to be the image
    set(hsurf,'CData', flipud(Z(:,:,i))');
    % set some extra properties
    set(hsurf,'EdgeColor','none');
    alpha(.5)
end

% Makeup
% set the viewing angle
view(g.comodtazimuth, g.comodtelevation)
axis tight; grid on;
h_ylabel = ylabel(g.ylabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
h_zlabel = zlabel(g.zlabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
h_xlabel = xlabel(g.xlabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');

set(haxes,'Color','None');

% Y axis labels
ylabel_val = freqs_amp(get(haxes,'YTick'));
for i = 1:length(ylabel_val)
    ylabel_string{i} = sprintf('%1.1f',ylabel_val(i));
end
set(haxes,'YTickLabel', ylabel_string,'FontSize',AXES_FONTSIZE_L);

% Z axis label
zlabel_val = freqs_phase(get(haxes,'ZTick'));
for i = 1:length(zlabel_val)
    zlabel_string{i} = sprintf('%1.1f',zlabel_val(i));
end
set(haxes,'ZTickLabel', zlabel_string);

% X axis label
set(haxes,'XTick', tplot);
for i = 1:length(tplot)
    xlabel_string{i} = sprintf('%1.2f',tplot(i));
end
set(haxes,'XTickLabel', xlabel_string);

% Title
title(g.title,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold');

% colorbar
hbar = colorbar;
set(get(hbar,'Label'),'String',g.cbartext);  
end
    