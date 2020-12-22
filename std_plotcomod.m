% std_plotcomod() - plot STUDY Comodulogram. 

% Usage:    
%   >> [STUDY] = std_erspplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY erspdata ersptimes erspfreqs pgroup pcond pinter] = ...
%                std_erspplot(STUDY, ALLEEG ...);
%
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
%   either 'channels' or 'cluster' inputs are also mandatory.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERSP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERSP of all subjects.
%   'noplot'   - ['on'|'off'] When 'on', only return output values. Default
%                is 'off'.
%
% Optional inputs:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERSPs of the requested clusters are plotted in the same figure, 
%                with ERSPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERSP for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERSP plus the
%                average cluster ERSP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Other optional inputs:
%   'plotmode'  - ['normal'|'condensed'|'none'] 'normal'  -> plot in a new figure; 
%                 'condensed' -> plot all curves in the current figure in a 
%                 condensed fashion. 'none' toggles off plotting {default: 'normal'}
%   'key','val' - All optional inputs to pop_specparams() are also accepted here
%                 to plot subset of time, statistics etc. The values used by default
%                 are the ones set using pop_specparams() and stored in the
%                 STUDY structure.
% Output:
%   STUDY      - the input STUDY set structure with the plotted cluster 
%                mean ERSPs added to allow quick replotting 
%   erspdata   - [cell] ERSP data for each condition, group and subjects.
%                size of cell array is [nconds x ngroups]. Size of each element
%                is [freqs x times x subjects] for data channels or 
%                [freqs x times x components] for component clusters. This 
%                array may be gicen as input  directly to the statcond() f
%                unction or std_stats() function to compute statistics.
%   ersptimes  - [array] ERSP time point latencies.
%   erspfreqs  - [array] ERSP point frequency values.
%   pgroup     - [array or cell] p-values group statistics. Output of the 
%                statcond() function.
%   pcond      - [array or cell] condition statistics. Output of the statcond() 
%                function.
%   pinter     - [array or cell] groups x conditions statistics. Output of
%                statcond() function.
%
% Example:
%        >> [STUDY] = std_erspplot(STUDY,ALLEEG, 'clusters', 'all', ...
%                                       'mode', 'together');
%           % Plot the mean ERSPs of all clusters in STUDY together 
%           % on the same figure. 
%
% Known limitations: when plotting multiple clusters, the output
%                    contains the last plotted cluster.
%
% See also: pop_clustedit(), pop_preclust(), eeg_createdata(), eeg_preclust(), pop_clustedit()
%
% Authors: Arnaud Delorme, CERCO, August, 2006

% Copyright (C) Arnaud Delorme, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [STUDY, allcomod, allfreqs1, allfreqs2, pgroup, pcond, pinter] = std_plotcomod(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_plotcomod;
    return;
end

STUDY  = pop_comodpacparams(STUDY, varargin{:});
params = STUDY.etc.pacplotopt.comodparam;
options = mystruct(varargin);
options = myrmfield( options, myfieldnames(params));

% get stat parameters
statstruct.etc = STUDY.etc; 
statstruct = pop_statparams(statstruct, varargin{:});

[ opt, moreparams ] = finputcheck( options, { ...
                               'pacmethod'      'string'        ''          '';
                               'design'      'integer' [] STUDY.currentdesign;
                               'caxis'       'real'    [] [];
                               'statmode'    'string'  [] ''; % deprecated
                               'channels'    'cell'    []              {};
                               'clusters'    'integer' []              [];
                               'datatype'    'string'  { 'pac','mipac' }  'pac';
                               'plottf'      'real'    []              [];
                               'mode'        'string'  []              ''; % for backward compatibility (now used for statistics)
                               'comps'       {'integer','string'}  []              []; % for backward compatibility
                               'plotsubjects' 'string' { 'on','off' }  'off';
                               'noplot'      'string'  { 'on','off' }  'off';
                               'plotmode'    'string'  { 'normal','condensed','none' }  'normal';
                               'subject'     'string'  []              '' }, ...
                               'std_plotcomod', 'ignore');

if ischar(opt), error(opt); end
if strcmpi(opt.noplot, 'on'), opt.plotmode = 'none'; end

allconditions  = {};
allgroups      = {};
condname       = '';
groupname      = '';
if length(STUDY.design(opt.design).variable) > 0, allconditions = STUDY.design(opt.design).variable(1).value; condname  = STUDY.design(opt.design).variable(1).label; end
if length(STUDY.design(opt.design).variable) > 1, allgroups     = STUDY.design(opt.design).variable(2).value; groupname = STUDY.design(opt.design).variable(2).label; end

stats = statstruct.etc.statistics;
stats.fieldtrip.channelneighbor = struct([]); % asumes one channel or 1 component
if isempty(STUDY.design(opt.design).variable)
    stats.paired = { };
else
    stats.paired = { STUDY.design(opt.design).variable(:).pairing };
end
if strcmpi(stats.singletrials, 'off') && ((~isempty(opt.subject) || ~isempty(opt.comps)))
    if strcmpi(stats.condstats, 'on') || strcmpi(stats.groupstats, 'on')
        stats.groupstats = 'off';
        stats.condstats  = 'off'; 
        disp('No statistics for single subject/component'); 
    end
end

if length(opt.comps) == 1
    stats.condstats = 'off'; stats.groupstats = 'off'; 
    disp('Statistics cannot be computed for single component');
end

alpha    = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.alpha, stats.fieldtrip.alpha);
mcorrect = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.mcorrect, stats.fieldtrip.mcorrect);
method   = fastif(strcmpi(stats.mode, 'eeglab'), stats.eeglab.method, ['Fieldtrip ' stats.fieldtrip.method ]);
plottfopt = { ...
   'ersplim',     opt.caxis, ...
   'threshold',   alpha, ...
   'effect',      stats.effect, ...
   'averagemode'  params.jointimemethod };

% Channels
if ~isempty(opt.channels)
        if length(opt.channels) > 1 && strcmpi(stats.singletrials, 'on')
        error('Cannot plot several channels on the same figure when using single trial statistics');
    end

%    [STUDY, allersp, alltimes, allfreqs, events, paramsersp] = std_readpacdata(STUDY, ALLEEG, 'channels', opt.channels, 'timerange', params.timerange,'freqrange', params.freqrange, 'subject', opt.subject, 'singletrials', stats.singletrials, 'design', opt.design, 'datatype', opt.datatype, 'subbaseline', params.subbaseline);
     [STUDY, allpacvals, alltimes, freq1, freq2, events, fileparams] = std_readpacdata(STUDY, ALLEEG, 'channels', opt.channels,         'timerange',params.timerange, ...
                                                                                                     'freqrange1', params.freqrange1,  'freqrange2', params.freqrange2,...
                                                                                                     'subject', opt.subject,           'singletrials',  stats.singletrials,...
                                                                                                     'design', opt.design,             'datatype',opt.datatype);
         if strcmpi(opt.datatype,'MIPAC') &&  strcmpi(stats.singletrials, 'off')
             % average MIPAC trials here
         end
         
          if strcmpi(opt.datatype,'PAC')
              for index = 1:length(allpacvals(:))
                  allpacvals{index} = squeeze(mean(allpacvals{index},3));
              end
         end
         
         

         [pcond, pgroup, pinter] = std_stat(allpacvals, stats);
         if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
                 (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1))
             pcond = {}; pgroup = {}; pinter = {};
             disp('No statistics possible for single subject STUDY');
         end % single subject STUDY
         
         locsOri = eeg_mergelocs(ALLEEG.chanlocs);
         locs = locsOri(std_chaninds(STUDY, opt.channels));


         if length(locs) > 1, opt.plottopo = 'on'; else opt.plottopo = 'off'; end
         if length(locs) == 1 && size(allpacvals{1},3) > 1
             % channels should be in 3rd dim; reshape data to put subjects in the 4th dim if number of channels is 1
             for index = 1:length(allpacvals(:))
                 allpacvals{index} = reshape(allpacvals{index}, size(allpacvals{index},1), size(allpacvals{index},2), 1, size(allpacvals{index},3));
             end
         end

         alltitles = std_figtitle('threshold', alpha, 'mcorrect', mcorrect, 'condstat', stats.condstats, 'cond2stat', stats.groupstats, ...
                                 'statistics', method, 'condnames', allconditions, 'cond2names', allgroups, ...
                                 'subject', opt.subject, 'datatype', upper(opt.datatype), 'plotmode', opt.plotmode, ...
                                 'effect', stats.effect, 'factor1', condname, 'factor2', groupname, 'chanlabels', opt.channels);
                             
         std_plotff(freq1, freq2, allpacvals, 'datatype', 'ersp', 'titles', alltitles, ...
             'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plottopo', 'off', 'plotmode', ...
             opt.plotmode, 'unitcolor', 'MI', 'chanlocs', locs, 'events', events);
else
    % Components(Clusters)
    if length(opt.clusters) > 1 && strcmpi(stats.singletrials, 'on')
        error('Cannot plot several components on the same figure when using single trial statistics');
    end
    
     nc = ceil(sqrt(length(opt.clusters)));
    nr = ceil(length(opt.clusters)/nc);
    comp_names = {};

    if length(opt.clusters) > 1 && ( strcmpi(stats.condstats, 'on') || strcmpi(stats.groupstats, 'on'))
        stats.condstats = 'off'; stats.groupstats = 'off';
    end
    
    for index = 1:length(opt.clusters)
        
        [STUDY, allpacvals, alltimes, freq1, freq2, events, fileparams] = std_readpacdata(STUDY, ALLEEG, 'clusters', opt.clusters(index), 'datatype', opt.datatype, ...
                                                                                                         'component', opt.comps,          'singletrials', stats.singletrials,...
                                                                                                         'timerange', params.timerange,   'freqrange1', params.freqrange1,...
                                                                                                         'freqrange2', params.freqrange2, 'design', opt.design);%,...
                                                                                                  %'concatenate', params.concatenate);
        
        if length(opt.clusters) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end
        
        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            comp_names = { STUDY.cluster(opt.clusters(index)).comps(opt.comps) };
            opt.subject = STUDY.datasetinfo(STUDY.cluster(opt.clusters(index)).sets(1,opt.comps)).subject;
        end
        
        % average single trials
        % ---------------------
        if strcmpi(opt.datatype, 'MIPAC')
            if  strcmpi(params.subbaseline, 'on')
                disp('Computing common baseline has changed since EEGLAB 14: averaging baselines is now');
                disp('performed before log-transformation of the baseline - in a similar way that baseline');
                disp('is averaged accross trials (log transformation is only performed at the end for display)');
                % see above for rational for baseline
                fileparams.singletrials = stats.singletrials;
                fileparams.commonbase   = params.subbaseline;
                allpacvals = newtimefbaseln(allpacvals, alltimes, fileparams);
            else
                fileparams.singletrials = stats.singletrials;
                allpacvals = cellfun(@(x)newtimefbaseln(x, alltimes, fileparams), allpacvals, 'uniformoutput', false);
            end
            % transform to log (except single trials)
            if strcmpi(stats.singletrials, 'off')
                if ~isfield(fileparams, 'scale') || strcmpi(fileparams.scale, 'log')
                    allpacvals = cellfun(@(x)10*log10(x), allpacvals, 'uniformoutput', false);
                end
            end
        end
        
                % statistics
        % ----------
        [pcond pgroup pinter] = std_stat(allpacvals, stats);
        if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
           (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1)), 
            pcond = {}; pgroup = {}; pinter = {}; 
            disp('No statistics possible for single subject STUDY');
        end % single subject STUDY                                
        
        if strcmpi(stats.singletrials, 'on')
            if strcmpi(opt.datatype, 'ersp')
                if ndims(allpacvals{1}) == 4, for ind = 1:length(allpacvals(:)), allpacvals{ind} = mean(allpacvals{ind},4); end; end
                if ndims(allpacvals{1}) == 3, for ind = 1:length(allpacvals(:)), allpacvals{ind} = mean(allpacvals{ind},3); end; end
                if ~isfield(fileparams, 'scale') || strcmpi(fileparams.scale, 'log')
                    allpacvals = cellfun(@(x)10*log10(x), allpacvals, 'uniformoutput', false);
                end
            elseif strcmpi(opt.datatype, 'itc')
                if ~isfield(params, 'itctype'), params.itctype = 'phasecoher'; end
                for iDat = 1:length(allpacvals(:))
                    allpacvals{iDat} = newtimefitc(allpacvals{iDat}, params.itctype);
                    allpacvals{iDat} = abs(allpacvals{iDat});
                end
            end
        end

        % plot specific component
        % -----------------------
        if index == length(opt.clusters), opt.legend = 'on'; end
        if ~strcmpi(opt.plotmode, 'none')
            alltitles = std_figtitle('threshold', alpha, 'mcorrect', mcorrect, 'condstat', stats.condstats, 'cond2stat', stats.groupstats, ...
                                     'statistics', method, 'condnames', allconditions, 'cond2names', allgroups, 'clustname', STUDY.cluster(opt.clusters(index)).name, 'compnames', comp_names, ...
                                     'subject', opt.subject, 'datatype', upper(opt.datatype), 'plotmode', opt.plotmode, ...
                                     'effect', stats.effect, 'factor1', condname, 'factor2', groupname);
            
            std_plotff(freq1, freq2, allpacvals, 'datatype', 'ersp', ...
                                           'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plotmode', ...
                                           opt.plotmode, 'titles', alltitles,'unitcolor', 'MI', 'chanlocs', ALLEEG(1).chanlocs, plottfopt{:});
        end
        
    end

end

% remove fields and ignore fields who are absent
% ----------------------------------------------
function s = myrmfield(s, f);

for index = 1:length(f)
    if isfield(s, f{index})
        s = rmfield(s, f{index});
    end
end
end
% convert to structure (but take into account cells)
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
end

% convert to structure (but take into account cells)
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
end
end