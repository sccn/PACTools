function hfig = eeg_plottfpac(EEG,varargin)

hfig = [];
if nargin < 1
    help eeg_plotfpact;
    return;
end

EEG  = pop_tfpacparams(EEG, varargin{:});
params = EEG.etc.pacplotopt.tfparam;
% Check validity of parameters

options = mystruct(varargin);
options = myrmfield( options, myfieldnames(params));

g  = finputcheck(  options, { ...
                   'plotindx'       'integer'       []          1;
                   'pacmethod'      'string'        ''          '';
                   'plotopt'        'cell'          {}          {}}, ...
                   'eeg_plotcomodt', 'ignore');

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
    disp('eeg_plottfpac() error: No PAC value has been computed for the input provided.');
    return;
end

pacdata   = pacstruct.(g.pacmethod).pacval;
freq1vals = pacstruct.params.freqs_phase;
freq2vals = pacstruct.params.freqs_amp;

if isfield( pacstruct.(g.pacmethod),'times')
    timevals  = pacstruct.(g.pacmethod).times;
else
    disp('eeg_plotfpac() message: Unable to generate plot. No time dimension found.');
    return;
end

if params.fixfreq == 1
    fixfreqval = params.freqval1;
else
    fixfreqval = params.freqval2;
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

% Trimming/fixing Phase frequency values
if params.fixfreq == 1
    if isempty(fixfreqval), disp('eeg_plotfpac() message: A value must be defined for the Phase frequency.'); return; end
    [~, indfreq1] = min(abs(freq1vals - fixfreqval));
    if isempty(indfreq1), error('eeg_plotfpac() message: Frequency value out of range'); end
    pacdata = pacdata(indfreq1,:,:,:);
else
    if ~isempty(params.freqrange)
        disp('Trimming Phase frequencies: Looking for the closest frequencies to the ones requested');
        indfreq1 = find(freq1vals > params.freqrange(1) & freq1vals < params.freqrange(2));
        freq1vals =  freq1vals(indfreq1);
        pacdata = pacdata(indfreq1,:,:,:);
        %     if flag_pval
        %         pvalmask = pvalmask(indfreq1,:,:,:);
        %     end
    end
end

% Trimming/fixing Amplitude frequency value

if params.fixfreq == 2
    if isempty(fixfreqval), disp('eeg_plotfpac() message: A value must be defined for the Phase frequency.'); return; end
    [~, indfreq2] = min(abs(freq2vals - fixfreqval));
    if isempty(indfreq2), error('eeg_plotfpac() message: Frequency value out of range'); end
    pacdata = pacdata(:,indfreq2,:,:);
else
    if ~isempty(params.freqrange)
        disp('Trimming Amplitude frequencies: Looking for the closest frequencies to the ones requested');
        indfreq2 = find(freq2vals > params.freqrange(1) & freq2vals < params.freqrange(2));
        freq2vals =  freq2vals(indfreq2);
        if isempty(indfreq2), error('eeg_plotfpac() message: Frequency value out of range'); end
        pacdata = pacdata(:,indfreq2,:,:);
        %     if flag_pval
        %         pvalmask = pvalmask(:,indfreq2,:,:);
        %     end
    end
end

% Collapsing Trials
if pacstruct.(g.pacmethod).dim == 3
    pacdata = mean(pacdata,3);    
end

% Trimming time dimesion
timevals = [];
if pacstruct.(g.pacmethod).dim == 2 || pacstruct.(g.pacmethod).dim == 3
    % Time values
    if  isfield(pacstruct.(g.pacmethod),'times')
        timevals = pacstruct.(g.pacmethod).times;
    end
        
    % Trimming time
    if ~isempty(params.timerange) && ~isempty(timevals)
        try
            timendx = (timevals > params.timerange(1) & timevals < params.timerange(2));
            timevals = timevals(timendx);
            if pacstruct.(g.pacmethod).dim == 2
                pacdata = pacdata(:,:,timendx);
            else
                pacdata = pacdata(:,:,:,timendx);
            end
        catch
            disp('Unable to trim time/latency values. Please check option ''timerange''');
            disp('Ignoring  option ''timerange'' ......');
        end
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

modes = {'Phase', 'Amplitude'} ;
if params.fixfreq == 1, 
     indxtitle = 2; 
     freqvals = freq2vals;
else indxtitle = 1;
     freqvals = freq1vals;
end

title  = [ modes{indxtitle} ' Vs Latency (f_{' modes{params.fixfreq} '} = ' num2str(fixfreqval) 'Hz)'];
TrialAxeTitle = [modes{indxtitle}  ' Frequency (Hz)'];

h = figure('Name', ['frequency-Latency PAC (' g.pacmethod ')'] ,'Units','Normalized'); hold on;
[~,~,~,~,axhndls] = erpimage(squeeze(pacdata)',[],timevals,title,0,0,'img_trialax_label',TrialAxeTitle,'cbar', 'on','cbar_title', 'Mod. Index','erp', 'on','yerplabel', 'Marginal Mod. Ind.', g.plotopt{:}) ;

set(get(axhndls{3},'Xlabel'),'String', 'Latency (ms)');
set(axhndls{1}, 'box', 'on'); set(axhndls{2}, 'box', 'on'); set(axhndls{3}, 'box', 'on');
set(axhndls{1},'YTickLabel',round(freqvals(get(axhndls{1},'YTick'))));
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