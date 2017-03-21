% eeg_visualize_pac() - visualize the phase amplitude coupling data
%
% Usage:
%   >> eeg_visualize_pac(EEG);
%   >> eeg_visualize_pac(EEG,'key1', 'val1', 'key2', val2' ...);
%
% Inputs:
%    EEG      - EEG structure or pac structure
%
%   Optional inputs
%       'time'              - Time point(s) to plot in seconds. Only one  
%                             time point is used for all the plots except
%                             plotcomodt. Default [1]                
%       'phasefreq'         - Frequency of the phase to plot. The plots will
%                             use the closest frequency available and the
%                             median frequency by default                           
%       'ampfreq'           - Frequency of the amplitude to plot. The plots will
%                             use the closest frequency available and the
%                             median frequency be default.
%       'plotcomod'         - Logical. Plot the comodulogram. 
%                             Default [1] if no other plot flags given,
%                             otherwise [0]
%       'plotcomodt'        - Logical. Plot the comodulogram at the given 
%                             time point(s). Several time points can be used.
%                             Default [0]
%       'plotampt'          - Logical. Plot the modulation index in time
%                             and in amplitudes at a given phase frequency.
%                             Default [0]
%       'plotphaset'        - Logical. Plot the modulation index in time
%                             and in phases at a given amplitude frequency.
%                             Default [0]
%       'plotsurrdist'      - Logical. Plot the surrogate distribution for 
%                             a given time point, amplitude frequency, 
%                             and phase frequency. Default [0]                          
%       'plotkl'            - Logical. Plot Kullback-Leibler modulation 
%                             index plot phase amplitude for a given time point,
%                             amplitude frequency, and phase frequency.                          
%                             Default [0]   
%       'plotmvl'           - Logical. Plot the distribution of the composites 
%                             of the mean vector length modulation index method
%                             for a given time point, amplitude frequency,
%                             and phase frequency. Default [0]
%       'arrowweight'       - Scalar used to multiply the ploted mean 
%                             vector length arrow. Default [0]
%       'nbinsmvl'          - Number of bins per phase bin and half the number 
%                             of phase bins to use in the composite distribution 
%                             of the mean vector length modulation index plot. 
%                             Default [36]
%       'normcomposite'     - Logical. Normalize the amplitude of the composites 
%                             values of the mean vector length modulation index. 
%                             Default [0].
%       'abspacval'         - Use the absolute value of the phase amplitude
%                             coupling value. Default [0]
%       'plotall'           - Logical. Plot all of the possible plots with
%                             the data and given parameters. Default [0]
%       'alphadata'         - Transparency to be used for to mask of
%                             significance in the comodulogram plot.
%                             Default [0.3]
%       'phasechanindx'     - When EEG structure is provided, will take
%                             this index to retreive and reconstruct the
%                             'pacstruct'
%       'ampchanindx'       - When EEG structure is provided, will take
%                             this index to retreive and reconstruct the
%                            'pacstruct'
% Outputs:
%         h                 - Cell array with the handles to the figures
%                             generated
%
% NOTE for developers: If changes in the inputs, update pop_pacplot as well
%                       
% Author: Joseph Heng and Ramon Martinez Cancino, EPFL, SCCN/INC, UCSD 2016
%
% See also:
%
% Copyright (C) 2002  Joseph Heng and Ramon Martinez Cancino, EPFL, UCSD, 
% INC, SCCN
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


function h = eeg_visualize_pac(EEG, varargin)
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

timeflag  = 0;
phaseflag = 0;
ampflag   = 0;
plotflag  = 0;

try g.time;                   timeflag    = 1; catch, g.time           = 1;              end;
try g.phasefreq;              phaseflag   = 1; catch,                                    end;  % Default set once pactruct is defined
try g.ampfreq;                ampflag     = 1; catch,                                    end;  % Default set once pactruct is defined
try g.plotcomod;              plotflag    = 1; catch, g.plotcomod      = 0;              end;
try g.plotcomodt;             plotflag    = 1; catch, g.plotcomodt     = 0;              end;
try g.plotampt;               plotflag    = 1; catch, g.plotampt       = 0;              end;
try g.plotphaset;             plotflag    = 1; catch, g.plotphaset     = 0;              end;
try g.plotsurrdist;           plotflag    = 1; catch, g.plotsurrdist   = 0;              end;
try g.plotkl;                 plotflag    = 1; catch, g.plotkl         = 0;              end;
try g.plotmvl;                plotflag    = 1; catch, g.plotmvl        = 0;              end;
try g.normcomposite;                           catch, g.normcomposite  = 0;              end;
try g.nbinsmvl;                                catch, g.nbinsmvl       = 36;             end;
try g.abspacval;                               catch, g.abspacval      = 0;              end;
try g.plotall;                                 catch, g.plotall        = 0;              end;
try g.alphadata;                               catch, g.alphadata      = 0.3;            end;  
try g.arrowweight;                             catch, g.arrowweight    = 1;              end;
try g.phasechanindx;                           catch, g.phasechanindx  = 1;             end;
try g.ampchanindx;                             catch, g.ampchanindx    = 1;             end;


if isempty(g.time),  timeflag = 0; end;
% Detect if EEG or pacstruct structure
if isfield(EEG,'etc')
    if isfield(EEG.etc,'eegpac') && ~isempty(EEG.etc.eegpac)

    % Reconstruct pacstruct from EEG.etc arguments using indices from
    % inputs for channels/components used fro the PAC cmputation
        if isempty(g.phasechanindx) || isempty(g.ampchanindx)
            error('eeg_visualize_pac() error: phasechanindx and ampchanindx must be provided');
        end       
        pacstruct = reconstruct_pacstruct(EEG, g.phasechanindx, g.ampchanindx);
        if isempty(pacstruct), return; end;
    else
        error('eeg_visualize_pac() error: Invalid or inexistent field eegpac');
    end
else
    pacstruct = EEG;
end

% Set default phase and/or amplitude frequency
if ~phaseflag
    g.phasefreq = median(pacstruct.freqs_phase);
end
if ~ampflag
    g.ampfreq = median(pacstruct.freqs_amp);
end
% Determine if the data is single strial or multiple triale
if numel(size(pacstruct.pacval)) == 2,
    s_trial = true;
else
    s_trial = false;
end

% Set default plot
if ~plotflag && g.plotcomod == 0,
    g.plotcomod = 1;   
end

if g.plotall,
    g.plotcomod = 1;       
    if s_trial,
        if phaseflag && ampflag,
            if ~isempty(pacstruct.surrogate_pac),
                g.plotsurrdist = 1;
            end
            switch pacstruct.method
                case 'mvlmi'
                    g.plotmvl = 1;
                case 'klmi'
                    g.plotkl = 1;
            end
        end
    else
        if phaseflag,
            g.plotampt = 1;
        end
        if ampflag,
            g.plotphaset = 1;
        end
        if timeflag,
            g.plotcomodt = 1;
        end
        if phaseflag && ampflag,
            if timeflag,            
                if ~isempty(pacstruct.surrogate_pac),
                    g.plotsurrdist = 1;
                end
                switch pacstruct.method
                    case 'mvlmi'
                        g.plotmvl = 1;
                    case 'klmi'
                        g.plotkl = 1;  
                end
            end
        end
    end             
end
    
    
% Check inputs
if g.plotcomodt|| g.plotphaset || g.plotampt
    if s_trial
        disp('eeg_visualize_pac() : Impossible to plot in time with single trial data');
        g.plotcomodt   = 0;
        g.plotphaset= 0;
        g.plotampt  = 0;
    end
end

if g.plotsurrdist,
    if isempty(pacstruct.surrogate_pac)
        disp('eeg_visualize_pac() : Impossible to plot surrogate phase amplitude coupling, no surrogate data available');
        g.plotsurrdist = 0;
    end
end

if g.plotmvl,
    if ~strcmp(pacstruct.method,'mvlmi')
        disp('eeg_visualize_pac() : Impossible to plot Mean Vector Length plot on data that has not used the Mean Vector Length method');
        g.plotmvl;
    end
end

if g.plotkl,
    if ~strcmp(pacstruct.method,'klmi')
        disp('eeg_visualize_pac() : Impossible to plot Kullback Leibler Phase Amplitude Plot on data that has not used the Kullback Leibler method');
        g.plotkl = 0;
    end
end

if ~s_trial && timeflag == 0 && ( g.plotcomodt || g.plotsurrdist || g.plotmvl || g.plotkl ),
    fprintf('eeg_visualize_pac() : No time given in input. Using first time point by default \n');
end
if ampflag == 0 && ( g.plotphaset || g.plotsurrdist || g.plotmvl || g.plotkl),
    fprintf('eeg_visualize_pac() : No frequency of amplitude given. Using median frequency of data by default \n');
end
if phaseflag == 0 && ( g.plotampt || g.plotsurrdist || g.plotmvl || g.plotkl),
    fprintf('eeg_visualize_pac() : No frequency of phase given. Using median frequency of data by default \n');
end

%% Configure
icadefs;

pacval = pacstruct.pacval;
if g.abspacval,
    pacval = abs(pacval);
end

if timeflag && ~s_trial,
    % Find closest time point in timesout
    [~, closest_time_idx] = min(abs(repmat(pacstruct.timesout,length(g.time),1)-repmat(g.time',1,length(pacstruct.timesout))),'',2);
else
    closest_time_idx = 1;
    if timeflag,
        disp('eeg_visualize_pac() : Impossible to use time with single trial data');
    end
end

closest_time = pacstruct.timesout(closest_time_idx);

% Find closest phase frequency in freqs_phase
tmp = abs(pacstruct.freqs_phase-g.phasefreq);
[~, closest_phasefreq_idx] = min(tmp);
closest_phasefreq = pacstruct.freqs_phase(closest_phasefreq_idx);

% Find closest amplitude frequency in freqs_amp
tmp = abs(pacstruct.freqs_amp-g.ampfreq);
[~, closest_ampfreq_idx] = min(tmp);
closest_ampfreq = pacstruct.freqs_amp(closest_ampfreq_idx);

plot_indx = 0;

%--------------------------------------------------------------------------
% PLOTTING starts here
%--------------------------------------------------------------------------

%% Comodulogram
if g.plotcomod
    plot_indx = plot_indx +1;
    h(plot_indx) = figure('Name', 'Modulation Index','Units','Normalized','Position', [3.4740e-01 2.4352e-01 3.4219e-01 5.4074e-01], 'Resize', 'off','Tag','comod_plot');
    h_axes1 = axes('parent',h(plot_indx) , 'Units', 'Normalized',...
                   'ActivePositionProperty', 'outerposition',...
                   'Color','None',...
                   'Visible', 'off');
    haxes = axes('Units', 'Normalized',...
                'Position',h_axes1.Position,...
                'ActivePositionProperty', 'outerposition');
    
    imagesc(pacstruct.freqs_phase,pacstruct.freqs_amp,mean(pacval,3)', 'parent', haxes);
    set(haxes,'YDir','normal');
    title(haxes, sprintf(' Modulation Index'));
    colorbar(haxes);
    xlabel(haxes,'Phase Frequency [Hz]');
    ylabel(haxes,'Amplitude Frequency [Hz]');
    set(get(h(plot_indx),'Children'),'Fontsize',AXES_FONTSIZE_L+5);

     % Plot significance mask 
    if ~isempty(pacstruct.signifmask) && s_trial, 
        haxes2 = axes('Position',haxes.Position,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
        haxes2.ActivePositionProperty = 'outerposition';
        linkaxes([haxes,haxes2]);
        h_trans = imagesc(pacstruct.freqs_phase,pacstruct.freqs_amp,pacstruct.signifmask', 'parent', haxes2);
        set(h_trans,'AlphaData',g.alphadata);
        set(haxes2,'Color','None','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','')
        colormap(haxes,'jet');
        colormap(haxes2,'gray');
    end
    box on; grid on;
end

%% Plot Comodulagram at timepoint

if g.plotcomodt
    if length(g.time) == 1
        plot_indx = plot_indx + 1;
        h(plot_indx) = figure('Name', 'Modulation Index','Units','Normalized','Position', [3.4740e-01 2.4352e-01 3.4219e-01 5.4074e-01], 'Resize', 'off','Tag','comodt_plot');
        h_axes1 = axes('Parent',h(plot_indx),'Units', 'Normalized',...
            'ActivePositionProperty', 'outerposition',...
            'Color','None',...
            'Visible', 'off');
        haxes = axes('Units', 'Normalized',...
            'Position',h_axes1.Position,...
            'ActivePositionProperty', 'outerposition');
        imagesc(pacstruct.freqs_phase,pacstruct.freqs_amp,pacval(:,:,closest_time_idx)', 'parent', haxes);
        
        set(haxes,'YDir','normal');
        title(haxes, sprintf('Modulation Index \n Time = %.3f s', closest_time));
        colorbar(haxes);
        xlabel(haxes,'Phase Frequency [Hz]');
        ylabel(haxes,'Amplitude Frequency [Hz]');
        set(get(h(plot_indx),'Children'),'Fontsize',AXES_FONTSIZE_L+5);
        
        % Plot significance mask
        if ~isempty(pacstruct.signifmask)
            haxes2 = axes('Position',haxes.Position,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
            haxes2.ActivePositionProperty = 'outerposition';
            linkaxes([haxes,haxes2]);
            h_trans = imagesc(pacstruct.freqs_phase,pacstruct.freqs_amp,pacstruct.pval(:,:,closest_time_idx)', 'parent', haxes2);
            set(h_trans,'AlphaData',g.alphadata);
            set(haxes2,'Color','None','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','')
            colormap(haxes,'jet');
            colormap(haxes2,'gray');
        end
        box on; grid on;
        %% Plot slices of modulation index in time
    else
        plot_indx = plot_indx + 1;
        h(plot_indx)     = figure('Name', 'Modulation Index in time','Units','Normalized','Position', [0.2349    0.3093    0.5547    0.2907],'Tag','comodt_plot');
        haxes = axes('Units', 'Normalized',...
            'Color','None',...
            'parent', h(plot_indx));
        
        % Plot
        indxval = closest_time_idx;
        Z       = pacval(:,:,indxval);
        tplot = pacstruct.timesout(closest_time_idx);
        [M,N,P] = size(Z);
        for i=1:P
            % Create a plane at x=i
            h(plot_indx) = surface(tplot(i)*ones(1,M),1:N,repmat([M:-1:1],N,1),repmat([M:-1:1],N,1));
            % set the color of the plane to be the image
            set(h(plot_indx),'CData', flipud(Z(:,:,i))');
            % set some extra properties
            set(h(plot_indx),'EdgeColor','none', 'FaceColor','interp')
            alpha(.5)
        end
        
        % Makeup
        % set the viewing angle
        view(3)
        axis tight; grid on;
        h_ylabel = ylabel('Amplitude Frequency (Hz)','FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
        h_zlabel = zlabel('Phase Frequency (Hz)','FontSize'    ,AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
        h_xlabel = xlabel('Time(sec)','FontSize'               ,AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
        
        set(get(h(plot_indx),'Parent'),'Color','None');
        
        % Amplitude axis labels
        ylabel_val = pacstruct.freqs_amp(get(get(h(plot_indx),'Parent'),'YTick'));
        for i = 1:length(ylabel_val)
            ylabel_string{i} = sprintf('%1.1f',ylabel_val(i));
        end
        set(get(h(plot_indx),'Parent'),'YTickLabel', ylabel_string,'FontSize',AXES_FONTSIZE_L);
        
        % Phase axis label
        zlabel_val = pacstruct.freqs_phase(get(get(h(plot_indx),'Parent'),'ZTick'));
        for i = 1:length(zlabel_val)
            zlabel_string{i} = sprintf('%1.1f',zlabel_val(i));
        end
        set(get(h(plot_indx),'Parent'),'ZTickLabel', zlabel_string);
        
        title('Modulation Index','FontSize',AXES_FONTSIZE_L,'FontWeight','bold');
        hbar = colorbar;
        set(get(hbar,'Label'),'String','Modulation Index');
    end
end

%% Plot modulation index in time
if g.plotphaset
    plot_indx = plot_indx + 1;
    h(plot_indx) = plot_pactimefreq(pacval,closest_ampfreq_idx,closest_ampfreq,pacstruct,1,'phaset_plot');
end
if g.plotampt  
    plot_indx = plot_indx + 1;
    h(plot_indx) = plot_pactimefreq(pacval,closest_phasefreq_idx,closest_phasefreq,pacstruct,0,'ampt_plot');
end

%% Plot surrogate distribution
if g.plotsurrdist
    surrogate_pactmp = squeeze(pacstruct.surrogate_pac(closest_phasefreq_idx, closest_ampfreq_idx,closest_time_idx, :));
    surrogate_pactmp = reshape(surrogate_pactmp,length(surrogate_pactmp),length(closest_time_idx));
    for itime = 1:length(closest_time_idx)
        surrogate_pac = surrogate_pactmp(:,itime);
        pacvaltmp = pacval(closest_phasefreq_idx, closest_ampfreq_idx, closest_time_idx(itime));
        p_value = pacstruct.pval(closest_phasefreq_idx, closest_ampfreq_idx, closest_time_idx(itime));
        if s_trial,
            figure_title = sprintf('Surrogate coupling value at %.1f Hz of amplitude frequency, %.1f Hz of phase frequency', closest_ampfreq, closest_phasefreq);
        else
            figure_title = sprintf('Surrogate coupling value at %.1f s, %.1f Hz of amplitude frequency, %.1f Hz of phase frequency', closest_time(itime), closest_ampfreq, closest_phasefreq);
        end
        h3 = figure('Name', figure_title, 'Units','Normalized', 'Position', [0.2456 0.2308 0.5244 0.5583],'Tag','surrdist_plot');
        haxes = axes('parent', h3);
        hhist = histfit(squeeze(surrogate_pac),max(1,ceil(length(surrogate_pac)/5)));
        set(hhist(1),'FaceColor',[0.4000    0.6980    1.0000]);
        set(hhist(2),'LineWidth',4);
        hold on;
        line([pacvaltmp pacvaltmp],haxes.YLim, 'Color', [1 1 1], 'parent', haxes);
        ylimtmp = get(haxes,'Ylim');
        xlimtmp = get(haxes,'Xlim');
        [x1fig y1fig] = axescoord2figurecoord(pacvaltmp, ylimtmp(2)/3, haxes);
        [x2fig y2fig] = axescoord2figurecoord(pacvaltmp, 0, haxes);
        harrow = annotation(h3, 'textarrow',[x1fig x2fig],[y1fig y2fig],'String','Observed PAC');
        set(harrow,'LineWidth',4,'FontSize',20,'HeadWidth',20,'HeadLength',20,'HeadStyle','plain');
        
        xlabel(haxes,'PAC Value');
        ylabel(haxes,'# surrogate PAC observations');
        if s_trial
            title(haxes, sprintf('Distribution of surrogate PAC \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz',...
                closest_phasefreq, closest_ampfreq));
        else
            title(haxes, sprintf('Distribution of surrogate PAC \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz  Time = %.1f s',...
                closest_phasefreq, closest_ampfreq, closest_time(itime)));
        end
        hlegend = legend(haxes,'Distribution of surrogates PAC values','Estimated normal distribution');
        set(hlegend,'box','off','FontSize',AXES_FONTSIZE_L+5);
        set(haxes,'FontSize',AXES_FONTSIZE_L+5);
        grid on;
        
        % Saving handles
        plot_indx = plot_indx + 1;
        h(plot_indx) = h3;
    end
end

%% Plot KLMI
if g.plotkl
    nbins = pacstruct.nbinskl(closest_phasefreq_idx, closest_ampfreq_idx, closest_time_idx);
    bin_average = pacstruct.bin_average{closest_phasefreq_idx, closest_ampfreq_idx, closest_time_idx};

    % Plot the bin_average
    h4 = figure('Name', 'Phase-Amplitude distribution','Units', 'Normalized','Position',[0.2995    0.3528    0.3698    0.3565],'Tag','klhist_plot');
    haxes = axes('parent',h4);
    bar(haxes, linspace(0,2*pi,nbins), bin_average);   
    set(haxes,'FontSize',AXES_FONTSIZE_L+5);
    xlabel(haxes,'Phase (rad)');
    ylabel(haxes,'Normalized Amplitude ');
    if s_trial
        title(haxes, sprintf('Phase-Amplitude distribution \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz', ...
                             closest_phasefreq, closest_ampfreq));
    else
        title(haxes, sprintf('Phase Amplitude Plot \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz \n Time point = %.3f s', ...
                             closest_phasefreq, closest_ampfreq, closest_time));            
    end
    xlim(haxes, [-2*pi/nbins, 2*pi+2*pi/nbins]);
    
    % Saving handles
    plot_indx = plot_indx + 1;
    h(plot_indx) = h4;
end

%% Plot MVLMI distribution
if g.plotmvl
    % Procesing for the plot
    n = g.nbinsmvl;
    compositestmp = squeeze(pacstruct.composites(closest_phasefreq_idx, closest_ampfreq_idx, closest_time_idx, :));
    compositestmp = reshape(compositestmp,length(compositestmp),length(closest_time_idx));
    for icomposites = 1: length(closest_time_idx)
        composites = compositestmp(:,icomposites);
        angle_composite   = wrapTo2Pi(angle(composites));
        amp_composite     = abs(composites);
        max_amp_composite = max(amp_composite);
        
        if g.normcomposite
            amp_composite = amp_composite/max_amp_composite;
            r_stepsize    = 1/n;
            r             = (0:n)'/n;
            mraw          = mean(amp_composite .* exp(1j*angle_composite));
        else
            r_stepsize    = max_amp_composite/n;
            r             = (0:r_stepsize:max_amp_composite)';
            mraw          = mean(composites);
        end
        
        theta = pi*(0:2*n)/n;
        theta_stepsize = pi/n;
        X = r*cos(theta);
        Y = r*sin(theta);
        
        C = zeros(length(r), length(theta));
        
        % Binning the composite values in the Real and Img part
        for i=1:length(r)
            for j=1:length(theta)
                C(i,j) = sum(angle_composite > wrapTo2Pi((j-1)*theta_stepsize) & ...
                    angle_composite < wrapTo2Pi(j*theta_stepsize)     & ...
                    amp_composite              > (i-1)* r_stepsize    & ...
                    amp_composite              <   i  * r_stepsize);
            end
        end
        
        % --- Plotting ---
        h6 = figure('Name', 'Distribution of Composite Values', 'Units', 'Normalized','Position',[0.3031 0.4241 0.3422 0.4806],'Tag','mvhist_plot');
        haxes = axes('parent', h6);
        pcolor(haxes, X,Y,C);
        hold on;
        hquiver = quiver(haxes, 0, 0, real(mraw), imag(mraw), g.arrowweight, 'r', 'linewidth', 3 ,'MaxHeadSize', 10);
        axis(haxes, 'equal', 'tight');
        shading(haxes,'flat');
        hcolorbar = colorbar(haxes);
        set(get(hcolorbar,'Label'),'String','# composite values')
        xlabel(haxes,'Real');
        ylabel(haxes,'Imaginary');
        if s_trial,
            title(haxes, sprintf('Distribution of Composite Values \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz',...
                closest_phasefreq, closest_ampfreq));
        else
            title(haxes, sprintf('Distribution of Composite Values \n Frequency_{phase} = %.1f Hz Frequency_{amplitude} = %.1f Hz \n Time = %.3f s',...
                closest_phasefreq, closest_ampfreq, closest_time(icomposites)));
        end
        set(haxes,'FontSize',AXES_FONTSIZE_L+5,'Color','None');
        hlegend = legend(hquiver,'Mean PAC vector');
        set(hlegend,'box', 'on','EdgeColor','None');
        set(get(hlegend,'BoxFace'),'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.5]))
        
        % Saving handles
        plot_indx = plot_indx + 1;
        h(plot_indx) = h6;
    end
end
end
%------------------------------------------------------------------------
function hfig = plot_pactimefreq(pacval,closest_freq_idx,closest_freq,pacstruct,ampflag,tagtext)
icadefs;

if ampflag
    ampphase = 'amplitude';
    mi = squeeze(pacval(:, closest_freq_idx, :));
    y_freqs = pacstruct.freqs_phase;
    y_label = 'Phase Frequencies [Hz]';
else
    ampphase = 'phase';
    mi = squeeze(pacval(closest_freq_idx, :, :));
    y_freqs = pacstruct.freqs_amp;
    y_label = 'Amplitude Frequencies [Hz]';
end

[~,~,ntimepoints] = size(pacval);
hfig = figure('Name', sprintf('Modulation Index for %.1f Hz of %s frequency', closest_freq,ampphase), 'Units','Normalized', 'Position', [0.2456 0.2308 0.5244 0.5583],'Tag',tagtext);
haxes = axes('parent',hfig,'Position',[0.1 0.21 8.500e-01 0.700]);
imagesc(pacstruct.timesout,y_freqs,mi);
colorbar(haxes);
set(haxes,'YDir','normal');
ylabel(haxes,y_label);
title(haxes, sprintf('Modulation Index \n Frequency_{%s} = %.1f Hz',ampphase,closest_freq));

% Axes 2
haxespos = get(haxes,'Position');
haxes = axes('parent',hfig,'Position',[haxespos(1)    0.09    haxespos(3)    0.12]);
plot(pacstruct.timesout,mean(mi), 'linewidth',3);
xlabel(haxes,'Time [s]');
ylabel(haxes,'Modulation Index');
set(get(hfig,'Children'),'FontSize',AXES_FONTSIZE_L+5);
axis(haxes, 'tight');
grid on;
end

function pacstruct = reconstruct_pacstruct(EEG, phasechanindx, ampchanindx)

    phaseidxstruct = find(EEG.etc.eegpac.phaseindx == phasechanindx);
    ampidxstruct   = find(EEG.etc.eegpac.ampindx   == ampchanindx);
    if all([~isempty(phaseidxstruct), ~isempty(ampidxstruct)])
    pacstruct                   = EEG.etc.eegpac.datapac{phaseidxstruct,ampidxstruct};
    pacstruct.method            = EEG.etc.eegpac.method;
    pacstruct.freqs_phase       = EEG.etc.eegpac.freqs_phase;
    pacstruct.freqs_amp         = EEG.etc.eegpac.freqs_amp;
    pacstruct.alpha             = EEG.etc.eegpac.alpha;
    pacstruct.ptspercent        = EEG.etc.eegpac.ptspercent;
    pacstruct.nboots            = EEG.etc.eegpac.nboots;
    pacstruct.srate             = EEG.etc.eegpac.srate;
    pacstruct.timesout          = EEG.etc.eegpac.timesout;
    pacstruct.options           = EEG.etc.eegpac.options;
%     pacstruct.freqs1            = EEG.etc.eegpac.freqs1;
%     pacstruct.freqs2            = EEG.etc.eegpac.freqs2;
    pacstruct.nbinskl           = EEG.etc.eegpac.nbinskl;
    else
        disp('eeg_visualize_pac() : Invalid data index');
        pacstruct = [];
        return;
    end
end