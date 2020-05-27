% std_checkifmeasures() - Check if the measure extension provided as...
%                         an input has been computed
% Usage:
%   >>  pac = std_checkifmeasures(STUDY, 'daterp');
% Inputs:
% STUDY     - EEGLAB STUDY stucture
% fileExt   - Extension of the STUDY measure to check e.g., 'daterp'
%
% Outputs:    
% AllMeasuresFlag  - [0|1] Binary variable informing of the presence [1] or
%                    not [0] of the measures being checked.
% MeasureFlag      - Binary vector with dimensions matching the numbers of
%                    subject in the STUDY informing of the presence [1] or
%                    not [0] of the measures being checked for each subject.
%                    Subject order is as in STUDY.datasetinfo.subject
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2019
%
% Copyright (C) 2020  Ramon Martinez-Cancino,INC, SCCN
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

function [AllMeasuresFlag, MeasureFlag] = std_checkifmeasures(STUDY,fileExt, varargin)

[opt moreopts] = finputcheck( varargin, {'design'        'integer' []             STUDY.currentdesign});
if ischar(opt), error(opt); end

% subjectList
allSubjects = { STUDY.datasetinfo.subject };
if isnan(opt.design)
    subjectList = STUDY.subject;
else
    subjectList = STUDY.design(opt.design).cases.value;
end
% get all sessions (same code as std_readdat)
% -------------------------------------------
allSessions = { STUDY.datasetinfo.session };
allSessions(cellfun(@isempty, allSessions)) = { 1 };
allSessions = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);

MeasureFlag = nan(1,length(subjectList));
for iSubj = 1:length(subjectList)
    datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
    fileName = getfilename({STUDY.datasetinfo(datInds).filepath}, STUDY.datasetinfo(datInds(1)).subject, { STUDY.datasetinfo(datInds).session }, ['.' fileExt], length(uniqueSessions) == 1);
    MeasureFlag(iSubj) = logical(exist(fileName,'file'));
end
AllMeasuresFlag = any(MeasureFlag);

%--------------------------------------------------------------------------|
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