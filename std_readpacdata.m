% std_reapacddata() - load pac measures for data channels or 
%                     for all components of a specified cluster.
%                     Called by plotting functions 
% Usage:
%         >> [STUDY, datavals, times] = ...
%                   std_readpacdata(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'design'    - [integer] read files from a specific STUDY design. Default
%                is empty (use current design in STUDY.currentdesign). Use
%                NaN to create a design with with all the data.
%  'channels'  - [cell] list of channels to import {default: none}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but 
%                the parent cluster (1) and any 'NotClust' clusters}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if 
%                available). Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'datatype'  - ['MIPAC'|'PAC'| select measure to load Default is 'PAC'.
%  'component' - [integer] select a specific component in a cluster.
%                This is the index of the component in the cluster not the
%                component number {default:all}
%  'timerange' - [min max] time range {default: whole measure range}
%
% Output:
%  STUDY    - updated studyset structure
%  pacvals  - [cell array] erp data (the cell array size is 
%             condition x groups)
%  times    - [float array] array of time values
%
% Example:
%  [pac times] = std_readpacdata(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(1).labels });
%
% Author: Ramon Martinez-Cancino, SCCN, UCSD
%         Arnaud Delorme, SCCN, UCSD
% Copyright (C) 2020  Ramon Martinez-Cancino, INC, SCCN
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

function [STUDY, datavals, xvals, yvals, zvals, events, params] = std_readpacdata(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readdata;
    return;
end

STUDY = pop_erspparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'design'        'integer'  []               STUDY.currentdesign;
    'channels'      'cell'     []               {};
    'clusters'      'integer'  []               [];
    'timerange'     'real'     []               [];
    'freqrange1'     'real'    []               [];
    'freqrange2'     'real'    []               [];
    'datatype'      'string'   { 'pac' 'mipac'} 'pac';
    'singletrials'  'string'   { 'on','off' }   'off';
    'componentpol'  'string'   { 'on','off' }   'on';
    'component'     'integer'  []               [];
    'subject'       'string'   []               '' }, ...
    'std_readdata', 'ignore');
if ischar(opt), error(opt); end

dtype = opt.datatype; % data type

% get the file extension
% ----------------------
tmpDataType = opt.datatype;

if ~isempty(opt.channels), fileExt = '.datpac';
else                       fileExt = '.icapac';
end

% list of subjects
% ----------------
allSubjects = { STUDY.datasetinfo.subject };
uniqueSubjects = unique(allSubjects);
STUDY.subject = uniqueSubjects;
if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end
if isempty(subjectList)
    if isnan(opt.design), subjectList = STUDY.subject;
    else subjectList = STUDY.design(opt.design).cases.value; 
    end
end

% options
% -------
opts = {};
if ~isempty(opt.timerange),  opts = [opts(:)'  {'timelimits'}, {opt.timerange} ];  end
if ~isempty(opt.freqrange1), opts = [opts(:)' {'freqlimits1'}, {opt.freqrange1} ]; end
if ~isempty(opt.freqrange2), opts = [opts(:)' {'freqlimits1'}, {opt.freqrange2} ]; end

opts = [ opts(:)' {'singletrials'} {opt.singletrials} ];
fprintf('Reading subjects'' data or looking up measure values in EEGLAB cache\n');

% get all sessions (same code as std_readdat)
% -------------------------------------------
allSessions = { STUDY.datasetinfo.session };
allSessions(cellfun(@isempty, allSessions)) = { 1 };
allSessions = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);

for iSubj = 1:length(subjectList)
    fprintf('.');
    
    % check cache
    bigstruct = [];
    if ~isempty(opt.channels), bigstruct.channel = opt.channels;
    else                       bigstruct.cluster = opt.clusters; % there can only be one cluster
    end
    bigstruct.datatype     = opt.datatype;
    bigstruct.singletrials = opt.singletrials;
    bigstruct.subject      = subjectList{iSubj};
    bigstruct.component    = opt.component;
    bigstruct.options      = opts;
    if isnan(opt.design)
         bigstruct.design.variable = struct([]);
    else bigstruct.design.variable = STUDY.design(opt.design).variable;
    end

    % find component indices
    % ----------------------
    if ~isempty(opt.clusters)
        datasetInds = strmatch(subjectList{iSubj}, { STUDY.datasetinfo.subject }, 'exact');
        compList    = [];
        polList     = [];
        if isempty(opt.component)
            for iDat = datasetInds(:)'
                indSet   = find(STUDY.cluster(opt.clusters).sets(1,:) == iDat); % each column contain info about the same subject so we many only consider the first row
                if ~isempty(indSet)
                    compList = [ compList STUDY.cluster(opt.clusters).comps(indSet) ];
                    if strcmpi(dtype, 'erp') && strcmpi(opt.componentpol, 'on')
                        polList  = [ polList  componentPol(indSet) ];
                    end
                end
            end
        else
            if ~isempty(intersect(datasetInds, STUDY.cluster(opt.clusters).sets(:,opt.component)))
                compList = [ compList STUDY.cluster(opt.clusters).comps(opt.component) ];
                if strcmpi(dtype, 'erp') && strcmpi(opt.componentpol, 'on')
                    polList  = [ polList  componentPol(opt.component) ];
                end
            end
        end
    end
    
    % read all channels/components at once
%     hashcode = gethashcode(std_serialize(bigstruct));
%     [STUDY.cache, tmpstruct] = eeg_cache(STUDY.cache, hashcode);
    tmpstruct = []; % RAMONcheck cache management for PAC
    if ~isempty(tmpstruct)
        dataTmp{iSubj}   = tmpstruct{1};
        xvals            = tmpstruct{2};
        yvals            = tmpstruct{3};
        zvals            = tmpstruct{4};
        eventsTmp{iSubj} = tmpstruct{5};
        params           = tmpstruct{6};
    else
        datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
        
        fileName = getfilename({STUDY.datasetinfo(datInds).filepath}, STUDY.datasetinfo(datInds(1)).subject, { STUDY.datasetinfo(datInds).session }, fileExt, length(uniqueSessions) == 1);
        if ~isempty(opt.channels)
             [dataTmp{iSubj}, params, xvals, yvals, zvals, eventsTmp{iSubj} ] = std_readpacfile( fileName, 'designvar', struct(bigstruct.design.variable), opts{:}, 'channels', opt.channels);
        else
            [dataTmp{iSubj}, params, xvals, yvals, zvals, eventsTmp{iSubj} ] = std_readpacfile( fileName, 'designvar', struct(bigstruct.design.variable), opts{:}, 'components', compList);
        end
        % Data manipulation removed here(--)
%         STUDY.cache = eeg_cache(STUDY.cache, hashcode, { dataTmp{iSubj} xvals yvals zvals eventsTmp{iSubj} params });
    end
end
fprintf('\n');

% if single trials, swap the last 2 dim (put channels before trials) %
% RAMON: Check when more than one channel
if strcmpi(opt.singletrials, 'on') && length(opt.channels) > 1
%     if ndims(dataTmp{1}{1}) == 4
%         for iCase = 1:length(dataTmp)
%             for iItem = 1:length(dataTmp{1}(:))
%                 dataTmp{iCase}{iItem} = permute(dataTmp{iCase}{iItem}, [1 3 2]);
%             end
%         end
%     else
  if ndims(dataTmp{1}{1}) == 5
        for iCase = 1:length(dataTmp)
            for iItem = 1:length(dataTmp{1}(:))
                dataTmp{iCase}{iItem} = permute(dataTmp{iCase}{iItem}, [1 2 3 5 4]);
            end
        end
  end
%     end
end

% store data for all subjects
if strcmp(opt.datatype, 'MIPAC')
    if length(opt.channels) > 1
        dim = 5; % 6
    else
        dim = 4;
    end
else
    if length(opt.channels) > 1
        dim = 5;
    else
        dim = 4; % ok for single channels no mipac
    end
end

events = {};
if ~isempty(opt.clusters)
    % Split ICA components from the same subjects need to be made 
    % as if coming from different subjects. E.g. 3 componets 2 S01 1S02 and
    % 2 cond -> will yield  dataTmp2 = {2×1 cell}    {2×1 cell}    {2×1 cell} 
    dataTmp2 = {};
    realDim  = dim;
    if strcmpi(opt.singletrials, 'on'), realDim = realDim+1; end
    for iDat1 = 1:length(dataTmp)
        compNumbers = cellfun(@(x)size(x, realDim), dataTmp{iDat1});
        if length(unique(compNumbers)) > 1
            error('Cannot handle conditions with different number of components');
        end
        
        for iComps = 1:compNumbers(1)
            dataTmp2{end+1} = [];
            for iDat2 = 1:length(dataTmp{iDat1}(:))
                % check dimensions of components
                if strcmpi(opt.singletrials, 'on') && strcmpi(tmpDataType, 'MIPAC'),    dataTmp2{end}{iDat2} = dataTmp{iDat1}{iDat2}(:,:,:,:,iComps);
                else                                                                    dataTmp2{end}{iDat2} = dataTmp{iDat1}{iDat2}(:,:,:,iComps);
                end
            end
            dataTmp2{end} = reshape(dataTmp2{end}, size(dataTmp{iDat1}));
        end
    end
    dataTmp = dataTmp2;
end
datavals = reorganizedata(dataTmp, dim); % here will pile up thee trials from all subjects in condition groups.

% reorganize data
% ---------------
function datavals = reorganizedata(dataTmp, dim)
    datavals = cell(size(dataTmp{1}));
        
    % copy data
    for iItem=1:length(dataTmp{1}(:)')
        numItems    = sum(cellfun(@(x)size(x{iItem},dim)*(size(x{iItem},1) > 1), dataTmp)); % the size > 1 allows to detect empty array which have a non-null last dim
        ind         = find(~cellfun(@(x)isempty(x{iItem}), dataTmp)); 
        if ~isempty(ind)
            ind = ind(1);
            switch dim
                case 2, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) numItems], 'single'); 
                case 3, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) numItems], 'single'); 
                case 4, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) size(dataTmp{ind}{iItem},3) numItems], 'single'); 
                case 5, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) size(dataTmp{ind}{iItem},3) size(dataTmp{ind}{iItem},4) numItems], 'single'); 
                case 6, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) size(dataTmp{ind}{iItem},3) size(dataTmp{ind}{iItem},4)  size(dataTmp{ind}{iItem},5) numItems], 'single'); 
            end
        end
    end
    for iItem=1:length(dataTmp{1}(:)')
        count = 1;
        for iCase = 1:length(dataTmp)
            if ~isempty(dataTmp{iCase}{iItem})
                numItems = size(dataTmp{iCase}{iItem},dim) * (size(dataTmp{iCase}{iItem},1) > 1); % the size > 1 allows to detect empty array which have a non-null last dim
                switch dim
                    case 2, datavals{iItem}(:,count:count+numItems-1) = dataTmp{iCase}{iItem}; 
                    case 3, datavals{iItem}(:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                    case 4, datavals{iItem}(:,:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                    case 5, datavals{iItem}(:,:,:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                    case 6, datavals{iItem}(:,:,:,:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                end
                count = count+numItems;
            end
        end
    end
    
% get file base name: filepath and sess are cell array (in case 2 files per subject)
% ----------------------------------------------------------------------------------
function filebase = getfilename(filepath, subj, sess, fileSuffix, onlyOneSession)
if onlyOneSession
    filebase = fullfile(filepath{1}, [ subj fileSuffix ] );
else
    if isempty(sess)
        sess = { '1' };
    end
    for iSess = 1:length(sess)
        if isnumeric(sess{iSess})
            sesStr   = [ '0' num2str(sess{iSess}) ];
        else
            sesStr   = [ '0' sess{iSess} ];
        end
        filebase{iSess} = fullfile(filepath{iSess}, [ subj '_ses-' sesStr(end-1:end) fileSuffix ] );
    end
end    