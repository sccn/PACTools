function hfig = eeg_plotcomod(EEG,varargin)

hfig = [];
if nargin < 1
    help eeg_plotcomod;
    return;
end

EEG  = pop_comodpacparams(EEG, varargin{:});
params = EEG.etc.pacplotopt.comodparam;
options = mystruct(varargin);
options = myrmfield( options, myfieldnames(params));

g  = finputcheck(  options, { ...
                   'plotsignif'     'integer'       [0 1]       0;
                   'plotindx'       'integer'       []          1;
                   'pacmethod'      'string'        ''          '';
                   'pval'           'integer'       ''          [];
                   'plotopt'        'cell'          {}          {}}, ...
                   'eeg_plotcomod', 'ignore');
                           
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
    disp('eeg_plotcomod() error: No PAC value has been computed for the input provided.');
    return;
end

pacdata   = pacstruct.(g.pacmethod).pacval;
freq1vals = pacstruct.params.freqs_phase;
freq2vals = pacstruct.params.freqs_amp;

% Pvals stuff
pvalmask = []; flag_pval =0;
if pacstruct.(g.pacmethod).dim ==1 && g.plotsignif  % Only for dim = 1 data 
    if isfield(pacstruct.(g.pacmethod),'signif') && ~isempty(pacstruct.(g.pacmethod).signif.pval)
        pvalmask = pacstruct.(g.pacmethod).signif.signifmask;
        flag_pval = 1;
    else
        disp('eeg_plotcomod message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
    end
end

if length(freq1vals)==1 || length(freq2vals)==1
    error('eeg_plotcomod() error: Insuficient number of frequencies for plotting');
end

% Trimming Phase frequency values
if ~isempty(params.freqrange1)
    if params.freqrange1(1) < min(freq1vals) || params.freqrange1(1) > max(freq1vals), error('eeg_plotcomod: Invalid phase frequency range'); end
    disp('Trimming Phase frequencies: Looking for the closest frequencies to the ones requested');
    indfreq1 = find(freq1vals > params.freqrange1(1) & freq1vals < params.freqrange1(2));
    freq1vals =  freq1vals(indfreq1);
    pacdata = pacdata(indfreq1,:,:,:);
    if flag_pval
        pvalmask = pvalmask(indfreq1,:,:,:);
    end
end

% Trimming Amplitude frequency values
if ~isempty(params.freqrange2)
    if params.freqrange2(1) < min(freq2vals) || params.freqrange2(2) > max(freq2vals), error('eeg_plotcomod: Invalid amplitude frequency range'); end
    disp('Trimming Amplitude frequencies: Looking for the closest frequencies to the ones requested');
    indfreq2 = find(freq2vals > params.freqrange2(1) & freq2vals < params.freqrange2(2));
    freq2vals =  freq2vals(indfreq2);
    pacdata = pacdata(:,indfreq2,:,:);
    if flag_pval
        pvalmask = pvalmask(:,indfreq2,:,:);
    end
end

% Collapsing Trials
if pacstruct.(g.pacmethod).dim == 3
    if strcmp(params.jointimemethod, 'average')
        pacdata = mean(pacdata,3);    
    elseif strcmp(params.jointimemethod, 'maxval')     
        pacdata = max(pacdata,[],3);   
    else
        disp('eeg_plotcomod() error: Invalid argument for property ''jointimemethod''. See ''pop_comodpacparams'' for resolution');
    end
end

% Trimming and collapsing time dimesion if present
timevals = [];
if pacstruct.(g.pacmethod).dim == 2 || pacstruct.(g.pacmethod).dim == 3
    % Time values
    if  isfield(pacstruct.(g.pacmethod),'times')
        timevals = pacstruct.(g.pacmethod).times;
    end
        
    % Trimming time
    if ~isempty(params.timerange) && ~isempty(timevals)
        if params.timerange(1) < min(timevals) || params.timerange(2) > max(timevals), error('eeg_plotcomod: Invalid time range'); end
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
    
    
    % Collapsing time
    if ~isempty(params.jointimemethod) && ~isempty(timevals)
        if strcmp(params.jointimemethod, 'average')
            if pacstruct.(g.pacmethod).dim == 2
                 pacdata = mean(pacdata,3);
            else pacdata = mean(pacdata,4);
            end
        elseif strcmp(params.jointimemethod, 'maxval')
            if pacstruct.(g.pacmethod).dim == 2
                 pacdata = max(pacdata,[],3);
            else pacdata = max(pacdata,[],4);
            end
        else
            disp('eeg_plotcomod() error: Invalid argument for property ''jointimemethod''. See ''pop_comodpacparams'' for resolution');
        end   
    end 
end

% Plotting significance
if flag_pval
    if ~isempty(pacstruct.(g.pacmethod).signif.signifmask)
        g.plotopt(end+1:end+2) = {'signifmask', pvalmask};
    else
        disp('eegplot pac() message: Significance mask overlay was requested but has not been computed. This input will be disregarded');
    end
end

if ~isempty(g.plotopt)
    postmp = ~cellfun(@isempty,strfind(g.plotopt(1:2:end),'title' ));
    if ~postmp
        g.plotopt{end+1} = 'title';
        g.plotopt{end+1} = ['Modulation Index ('  g.pacmethod ')'];
    end 
else
    g.plotopt = {'title', ['Modulation Index ('  g.pacmethod ')']}; 
end


hfig = comodulogram(freq1vals,freq2vals,pacdata, g.plotopt{:});

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