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