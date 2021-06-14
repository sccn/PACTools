function h = comodulogram(freqs_phase,freqs_amp,pacval,varargin)
icadefs;
h = [];
if nargin < 3
    help comodulogram;
    return;
end

g = finputcheck( varargin, { ...
                            'signifmask'     'integer'       [ ]                     [];
                            'alphadata'      'integer'       []                      0.5;
                            'title'          'string'        ''                      'Modulation Index';
                            'shading'        'string'        {'flat', 'faceted'}     'flat';
                            'scale'          'string'        {'linear' 'log'}        'log'});

h = figure('Name', g.title,'Units','Normalized','Tag','comod_plot');
haxes = axes('parent',h,'Units', 'Normalized','ActivePositionProperty', 'outerposition');
pcolor(freqs_phase, freqs_amp, pacval');
% shading(haxes,g.shading); % shading data
colormap(haxes,'parula');  % Colormap

% Scale
if strcmp(g.scale, 'log')
    set(haxes, 'XScale', 'log');
    set(haxes, 'YScale', 'log');
else
    set(haxes, 'XScale', 'linear');
    set(haxes, 'YScale', 'linear');
end

title(haxes, sprintf(g.title)); % Title
colorbar(haxes); % colorbar

% Axis labels
xlabel(haxes,'Phase Frequency [Hz]');
ylabel(haxes,'Amplitude Frequency [Hz]');
set(get(h,'Children'),'Fontsize',AXES_FONTSIZE_L+5);

% Plot significance mask
if ~isempty(g.signifmask)
    haxes2 = axes('Position',haxes.Position,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
    haxes2.ActivePositionProperty = 'outerposition';
    haxes2.Visible = 'off';
    linkaxes([haxes,haxes2]);
    h_trans = pcolor(freqs_phase, freqs_amp,int8(g.signifmask)', 'parent', haxes2);
    
    if strcmp(g.scale, 'log')
        set(haxes2, 'XScale', 'log');
        set(haxes2, 'YScale', 'log');
    else
        set(haxes2, 'XScale', 'linear');
        set(haxes2, 'YScale', 'linear');
    end
    set(h_trans,'FaceAlpha',g.alphadata);
    set(haxes2,'Color','None','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','')
    colormap(haxes,'parula');
    colormap(haxes2,'gray');
end
box on; grid on;
end