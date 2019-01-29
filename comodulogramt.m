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
try g.times;                           catch, g.times = [];                           end;
try g.title;                           catch, g.title = 'Modulation Index';           end;
try g.zlabel;                          catch, g.zlabel = 'Phase Frequency (Hz)';      end;
try g.ylabel;                          catch, g.ylabel = 'Amplitude Frequency (Hz)';  end;
try g.xlabel;                          catch, g.xlabel = 'Time(sec)';                 end;
try g.cbartext;                        catch, g.cbartext = 'Modulation Index';        end; 
try g.comodtazimuth;                   catch, g.comodtazimuth = -10;                  end;
try g.comodtelevation;                 catch, g.comodtelevation = 23;                 end;

% Determining nearest time
if isempty(g.times)
    times = min(timevect):((max(timevect)-min(timevect))/(g.ntimepoints-1)):max(timevect);
else
    times = g.times;
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