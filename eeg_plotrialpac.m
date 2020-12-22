function hfig = eeg_plotrialpac(EEG,varargin)

hfig = [];
if nargin < 1
    help eeg_plottrialpac;
    return;
end

EEG  = pop_trialspacparams(EEG, varargin{:});
params = EEG.etc.pacplotopt.trialparam;
% Check validity of parameters

options = mystruct(varargin);
options = myrmfield( options, myfieldnames(params));

g  = finputcheck(  options, { ...
                   'plotindx'       'integer'       []          1;
                   'pacmethod'      'string'        ''          '';
                   'plotopt'        'cell'          {}          {}}, ...
                   'eeg_plottrialpac', 'ignore');

pacstruct = EEG.etc.eegpac(g.plotindx);

% if not pacmethod provided use the first available
if isempty(g.pacmethod)
    fieldnameftmp = fieldnames(pacstruct);
    tmpval = find(structfun(@(x) ~isempty(x),EEG.etc.eegpac(g.plotindx)));
    tmpval2 = find(tmpval>5);
    g.pacmethod = fieldnameftmp{tmpval(tmpval2)};
end

% pacval may be empty
if isempty(pacstruct.(g.pacmethod).pacval)
    disp('eeg_plottrialpac() error: No PAC value has been computed for the input provided.');
    return;
end

% Only for dim == 3 data
if pacstruct.(g.pacmethod).dim ~= 3
    disp('eeg_plottrialpac() error: Invalid data for eeg_plottrialpac');
    return;
end

% Getting data
pacdata   = pacstruct.(g.pacmethod).pacval;
freq1vals = pacstruct.params.freqs_phase;
freq2vals = pacstruct.params.freqs_amp;

if isfield(pacstruct.(g.pacmethod),'times') && ~isempty(pacstruct.(g.pacmethod).times)
    timevals  = pacstruct.(g.pacmethod).times;
else
    disp('eeg_plotfpac() message: Unable to generate plot. No time dimension found.');
    return;
end

% % Pvals stuff
% pvalmask = []; flag_pval =0;
% if pacstruct.(g.pacmethod).dim ==1 && g.plotsignif  % Only for dim = 1 data 
%     if isfield(pacstruct.(g.pacmethod),'signif') && ~isempty(pacstruct.(g.pacmethod).signif.pval)
%         pvalmask = pacstruct.(g.pacmethod).signif.signifmask;
%         flag_pval = 1;
%     else
%         disp('eeg_plotcomod message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
%     end
% end

if length(freq1vals)==1 || length(freq2vals)==1
    error('eeg_plotfpac() error: Insuficient number of frequencies for plotting');
end

% Fixing Phase frequency values
if isempty(params.freqval1)
    nfreqstmp = numel(EEG.etc.eegpac(1).params.freqs_phase);
    params.freqval1 = EEG.etc.eegpac(1).params.freqs_phase(floor(nfreqstmp/2));
end

[~, indfreq1] = min(abs(freq1vals - params.freqval1));
if isempty(indfreq1), error('eeg_plotfpac() message: Frequency value out of range'); end
pacdata = pacdata(indfreq1,:,:,:);

% Trimming/fixing Amplitude frequency value
if isempty(params.freqval2)
    nfreqstmp = numel(EEG.etc.eegpac(1).params.freqs_amp);
    params.freqval2 = EEG.etc.eegpac(1).params.freqs_amp(floor(nfreqstmp/2));
end

[~, indfreq2] = min(abs(freq2vals - params.freqval2));
if isempty(indfreq2), error('eeg_plotfpac() message: Frequency value out of range'); end
pacdata = pacdata(:,indfreq2,:,:);

% Trimming time dimension
if ~isempty(params.timerange)
    try
        timendx = (timevals > params.timerange(1) & timevals < params.timerange(2));
        timevals = timevals(timendx);
        pacdata = pacdata(:,:,:,timendx);
    catch
        disp('Unable to trim time/latency values. Please check option ''timerange''');
        disp('Ignoring  option ''timerange'' ......');
    end
end

% Trimming trial dimesion
if ~isempty(params.trialrange)
    try  
        pacdata = pacdata(:,:,params.trialrange,:);
    catch
        disp('Unable to trim trials. Please check option ''trialrange''');
        disp('Ignoring  option ''trialrange'' ......');
    end
end

% % Plotting significance
% if flag_pval
%     if ~isempty(pacstruct.(g.pacmethod).signif.signifmask)
%         g.plotopt(end+1:end+2) = {'signifmask', pvalmask};
%     else
%         disp('eegplot pac() message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
%     end
% end
%---------------------------

title  = ['Trials Vs Latency (f_{Phase} = ' num2str(params.freqval1,'%.1f') 'Hz, f_{Amp} = ' num2str(params.freqval2,'%.1f') 'Hz )'];
plotopttmp = g.plotopt;
g.plotopt(end+1:end+length(plotopttmp)) = plotopttmp;
h = figure('Name', ['Trial-Time PAC (' g.pacmethod ')'] ,'Units','Normalized'); hold on;
[trash,trash,trash,trash,axhndls] = erpimage(squeeze(pacdata)',[],timevals,title,0,0,'img_trialax_label','Trials',...
                                                                            'cbar', 'on',...
                                                                            'cbar_title', 'Mod. Index',...
                                                                            'erp', 'on',...
                                                                            'yerplabel', 'Marginal Mod. Ind.', g.plotopt{:}) ;

set(get(axhndls{3},'Xlabel'),'String', 'Latency (ms)');
set(axhndls{1}, 'box', 'on'); set(axhndls{2}, 'box', 'on'); set(axhndls{3}, 'box', 'on');
set(axhndls{1},'Fontsize', 20);
set(axhndls{2},'Fontsize', 15);
set(axhndls{3},'Fontsize', 20);
htmp = get(axhndls{3},'Children');
set(htmp(1),'Fontsize', 15);

%--------------------------------------------------------------------------
% remove fields and ignore fields who are absent -------> std_erspplot
% ----------------------------------------------
function s = myrmfield(s, f);

for index = 1:length(f)
    if isfield(s, f{index})
        s = rmfield(s, f{index});
    end
end

% convert to structure (but take into account cells) -------> std_erspplot
% --------------------------------------------------
function s = mystruct(v);

for index=1:length(v)
    if iscell(v{index})
        v{index} = { v{index} };
    end
end
try
    s = struct(v{:});
catch, error('Parameter error'); end

% convert to structure (but take into account cells) -------> std_erspplot
% --------------------------------------------------
function s = myfieldnames(v);

s = fieldnames(v);
if isfield(v, 'eeglab')
    s2 = fieldnames(v.eeglab);
    s = { s{:} s2{:} };
end
if isfield(v, 'fieldtrip')
    s3 = fieldnames(v.fieldtrip);
    for index=1:length(s3)
        s3{index} = [ 'fieldtrip' s3{index} ];
    end
    s = { s{:} s3{:} };
end