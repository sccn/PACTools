function h = comodulogramt(freqs_phase,freqs_amp,timevect, pacval, varargin)

icadefs;
h = [];
if nargin < 4
    help comodulogramt;
    return;
end

DefaultTitle     = 'Modulation Index';
ZLabelDefault    = 'Amplitude Frequency (Hz)';
YLabelDefault    = 'Phase Frequency (Hz)';
XLabelDefault    = 'Latency(msec)';
CBarLabelDefault = 'Modulation Index';
g = finputcheck( varargin, { ...
                            'signifmask'      'integer'       []                         [];
                            'npoints'         'integer'       [2 length(timevect) ]      20;
                            'times'           'integer'       []                         [];
                            'title'           'string'        ''                         DefaultTitle;
                            'zlabel'          'string'        ''                         ZLabelDefault;
                            'ylabel'          'string'        ''                         YLabelDefault;
                            'xlabel'          'string'        ''                         XLabelDefault;
                            'cbartext'        'string'        ''                         CBarLabelDefault;
                            'comodtazimuth'   'integer'       []                         -10;
                            'comodtelevation' 'integer'       []                         23;
                            'facealpha'       'integer'       []                         0.5;
                            'scale'           'string'        {'linear' 'log'}           'log'});
                        
flag_pval = 0;
if ~isempty(g.signifmask) && ~isequal(size(g.signifmask), size(pacval))
    error('comodulogramt : PAC values and significnace mask dimensions diverge');
elseif ~isempty(g.signifmask) && isequal(size(g.signifmask), size(pacval))
    flag_pval = 1    
end
% Determining nearest time
if isempty(g.times)
    times = min(timevect):((max(timevect)-min(timevect))/(g.npoints-1)):max(timevect);
else
    times = g.times;
end
[~, timeidx] = min(abs(repmat(timevect,length(times),1)-repmat(times',1,length(timevect))),'',2);
tplot = timevect(timeidx);   

%% Plotting start here

h     = figure('Name', g.title,'Units','Normalized','Position', [0.2349    0.3093    0.5547    0.2907],'Tag','comodt_plot');
haxes = axes('Units', 'Normalized','Color','None','parent', h, 'YScale', g.scale, 'ZScale', g.scale);

% Plot
Z = permute(pacval(:,:,timeidx), [2,1,3]);
[M,N,P] = size(Z);
for i=1:P
    % Create a plane at x=i
    hsurf = surface(tplot(i)*ones(1,M),freqs_phase ,repmat(flip(freqs_amp),N,1),repmat([M:-1:1],N,1));
    
    if flag_pval
        set(hsurf,'CData',real(flipud(g.signifmask(:,:,i))').*flipud(Z(:,:,i))', 'EdgeColor','none');
    else
    % set the color of the plane to be the image
    set(hsurf,'CData', flipud(Z(:,:,i))', 'EdgeColor','none');
    end
alpha(g.facealpha)
end

% Makeup
% set the viewing angle
view(g.comodtazimuth, g.comodtelevation)
axis tight; grid on;
h_ylabel = ylabel(g.ylabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
h_zlabel = zlabel(g.zlabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');
h_xlabel = xlabel(g.xlabel,'FontSize',AXES_FONTSIZE_L,'FontWeight','bold','Units','Normalized');

set(haxes,'Color','None');

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
    