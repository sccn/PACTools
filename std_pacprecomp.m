% std_pacprecomp() - Compute PAC measures in the EEGLAB STUDY
%
% Usage:  
%   >> [STUDY, ALLEEG] = std_pac(STUDY, ALLEEGY, 'key', 'val', ...);
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%
% Optional inputs:
%
%
% Other optional inputs:
%   This function will take any of the newtimef() optional inputs (for instance
%   to compute log-space frequencies)...
%
% Outputs:
%

function [STUDY ALLEEG] = std_pacprecomp(STUDY, ALLEEG,freqs1, freqs2, chanorcomp,dataindx, varargin)

if nargin < 1
    help std_pacprecomp;
    return;
end
pacmethods_list = {'plv','mvlmi','klmi','glm','plv', 'instmipac', 'ermipac'} ;

[g pacargs] = finputcheck(varargin, { ...
                        'verbose'        'integer'     [0 1]               [1];
                        'savefile'       'string'      {'on', 'off'}       'on';
                        'datatype'       'integer'     [1 2]               [];
                        'trialindices'    { 'integer','cell' }     []      [];
                        'design'         'integer'     []                  [1];
                        'method'         'string'      pacmethods_list     'ermipac'; 
                        'timerange'      'real'        []                  [];  
                        'nsgopt'         'cell'        {}                  {'nsgflag', 0};
                        'compflag'       'string'      {'local', 'nsg'}    'local';
                        'pacparams'      'cell'        {}                  {};
                        'recompute'      'string'      {'on','off'}        'off';}, 'std_pac', 'ignore');
if ischar(g), error(g); end


%-------------------- Local computation --------------------
if strcmpi(g.compflag, 'local')
    dataprefix = {'dat', 'ica'};
    fieldprefix = {'chan', 'comp'};
    prefix = dataprefix{chanorcomp}; filefieldprefix = fieldprefix{chanorcomp};
    
    datatype = fastif(strcmp(g.method,'ermipac'),'MIPAC','PAC');
    
    % Check channels across subjects
    % here
    if isnumeric(dataindx)
        datindxint = dataindx;
        dataindx = {ALLEEG(1).chanlocs(dataindx).labels};
    else
        datindxint = find(contains({ALLEEG(1).chanlocs.labels},dataindx));
    end
    
    % Check for TF decomposition in STUDY
    AllMeasuresFlag = std_checkifmeasures(STUDY,[dataprefix{chanorcomp} 'timef']);
    if ~AllMeasuresFlag
        disp('std_pacprecomp: ERSP measure should be computed before use std_pacprecomp');
        exit;
    end
    
    % Check consistency of number of data points
    allPnts = [ALLEEG(:).pnts];
    if iscell(allPnts)
        allPnts = [ allPnts{:} ];
    end
    if length(unique(allPnts)) > 1
        error([ 'Cannot compute PAC in the current STUDY set because datasets' 10 'do not have the same number of data points' ])
    end
    
    % Extract PAC parameters
    if ~isempty(g.pacparams)
        optpac = cell2struct(g.pacparams(2:2:end),g.pacparams(1:2:end),2);
    end
    
    % Getting subjects and sessions
    allSubjects = { STUDY.datasetinfo.subject };
    allSessions = { STUDY.datasetinfo.session };
    uniqueSubjects = unique(allSubjects);
    allSessions(cellfun(@isempty, allSessions)) = { 1 };
    allSessions = cellfun(@num2str, allSessions, 'uniformoutput', false);
    uniqueSessions = unique(allSessions);
    
    if strcmpi(prefix, 'ica')
        clsindx = find(ismember({STUDY.cluster.name},dataindx));
        if isempty(clsindx) || all(clsindx ==0), error('Invalid input chanindx'), end
        subjICs = cell(1,length(uniqueSubjects));
        subjICs_indx = subjICs;
        subjIC_clust = subjICs;
        for icls = 1:length(clsindx)
            for iSubj = 1:length(uniqueSubjects)
                inds = strmatch( uniqueSubjects{iSubj}, allSubjects, 'exact');
                [~,J] = ind2sub(size(STUDY.cluster(clsindx(icls)).sets), find( STUDY.cluster(clsindx(icls)).sets==inds(1)));
                if ~isempty(J)
                    subjICs_indx{iSubj} = [subjICs_indx{iSubj} J'];
                    subjICs{iSubj} = [ subjICs{iSubj} STUDY.cluster(clsindx(icls)).comps(J)];
                    subjIC_clust{iSubj}  = [subjIC_clust{iSubj}  repmat(clsindx(icls),1,length(J))];
                end
            end
        end
        % Sorting components
        for iSubj = 1:length(uniqueSubjects)
            [subjICs{iSubj}, sortindx] = sort(subjICs{iSubj});
            subjICs_indx{iSubj} = subjICs_indx{iSubj}(sortindx);
            subjIC_clust{iSubj} = subjIC_clust{iSubj}(sortindx);
        end
        
    end
    
    for iSubj = 1:length(uniqueSubjects)
        for iSess = 1:length(uniqueSessions)
            
            inds = strmatch( uniqueSubjects{iSubj}, allSubjects, 'exact');
            filepath = STUDY.datasetinfo(inds(1)).filepath;
            trialinfo = std_combtrialinfo(STUDY.datasetinfo, inds);
            filebase = getfilename({filepath}, uniqueSubjects{iSubj}, uniqueSessions{iSess}, ['.' prefix 'pac'], length(uniqueSessions) == 1);
            
             if isempty(g.trialindices), g.trialindices = cell(length(ALLEEG(inds))); end
             if ~iscell(g.trialindices), g.trialindices = { g.trialindices }; end
            
            if strcmpi(prefix, 'dat')
                stream_indx = dataindx;
                % Checking consistency across datasets
                chan_indx = eeg_chaninds(ALLEEG(inds(1)), dataindx, 0);
                for datindx = 2:length(ALLEEG(inds))
                    if ~isequal(eeg_chaninds(ALLEEG(inds(datindx)), dataindx, 0), chan_indx)
                        error([ 'Channel information must be consistant when ' 10 'several datasets are merged for a specific design' ]);
                    end
                end
            else
                stream_indx = subjICs_indx{iSubj};
                cls_indx = subjIC_clust{iSubj};
            end
            
            
            for ichan = 1:length(stream_indx)
                if strcmpi(prefix, 'dat')
                    % Read TF decomposition for numchandinx(iChanIndx,1) and numchandinx(iChanIndx,2)
                    % Reading TF decomposition
                    % channels must be a cell
                    [~, alltf, alltimes, allfreqs, ~, paramstf] = std_readdata(STUDY, ALLEEG, 'channels', stream_indx(ichan),...
                                                                                'timerange', g.timerange, ...
                                                                                'freqrange', [min(freqs1) max(freqs2)],...
                                                                                'subject', uniqueSubjects{iSubj},...
                                                                                'singletrials', 'on',...
                                                                                'design', NaN,...
                                                                                'datatype', 'itc',...
                                                                                'subbaseline', 'off');
                else
                    % Channels must be a real
                    [~, alltf, alltimes, allfreqs, ~, paramstf] = std_readdata(STUDY, ALLEEG, 'component', stream_indx(ichan),...
                                                                                'clusters', cls_indx(ichan),...
                                                                                'timerange', g.timerange, ...
                                                                                'freqrange', [min(freqs1) max(freqs2)],...
                                                                                'subject', uniqueSubjects{iSubj},...
                                                                                'singletrials', 'on',...
                                                                                'design', NaN,...
                                                                                'datatype', 'itc',...
                                                                                'subbaseline', 'off');
                end
                
                % Select frequencies and narrow data to it
                if g.verbose
                    disp('Selecting closest frequencies to the ones provided.')
                end
                
                % Phase
                [~, indxfrqs1] = min(abs(repmat(allfreqs,length(freqs1),1)-repmat(freqs1,length(allfreqs),1)'),[],2);
                realfreqs1 =  allfreqs(indxfrqs1);
                for ifreqs =1:length(alltf)
                    allersp_1{ifreqs} =  alltf{ifreqs}(indxfrqs1,:,:);
                end
                
                if length(unique(indxfrqs1))~=length(freqs1)
                    sprintf(['std_pacprecom: Frequency values are not consistent with the computed ITC measure\n', ...
                        'This may happen when the frequency resolution of the input vector is higher than the\n',...
                        'one already computed. Consider, Modifying your frequency input vector or recompute the ITC\n'])
                    return;
                end
                alltf_1.alltf = allersp_1{1};
                alltf_1.freqs = realfreqs1;
                alltf_1.timesout = alltimes;
                
                % Amplitude
                [~, indxfrqs2] = min(abs(repmat(allfreqs,length(freqs2),1)-repmat(freqs2,length(allfreqs),1)'),[],2);
                realfreqs2 =  allfreqs(indxfrqs2);
                for ifreqs =1:length(alltf)
                    allersp_2{ifreqs} =  alltf{ifreqs}(indxfrqs2,:,:);
                end
                
                if length(unique(indxfrqs2))~=length(freqs2)
                    sprintf(['std_pacprecom: Frequency values are not consistent with the computed ITC measure\n', ...
                        'This may happen when the frequency resolution of the input vector is higher than the\n',...
                        'one already computed. Consider, Modifying your frequency input vector or recompute the ITC\n'])
                    return;
                end
                alltf_2.alltf = allersp_2{1};
                alltf_2.freqs = realfreqs2;
                alltf_2.timesout = alltimes;
                
                % Compute PAC
                [pacval, timesout, freqs1, freqs2, alltfX, alltfY,~, pacstruct] =...
                    eeg_pac([], [], ALLEEG(1).srate, 'alltfXstr',alltf_1,...
                    'alltfYstr',alltf_2,...
                    'method',g.method,...
                    'tlimits', alltimes, g.pacparams{:});
                if ndims(pacval)==4
                    pacval = permute(pacval,[1 2 4 3]);
                end
   
                allTrialsPAC{ichan}   = pacval;
            end
        end
        
        % Store PAC in structure
        all_pac = [];
        for ichan = 1:size(stream_indx,2)
            if strcmp(prefix,'dat')
%                 fieldname = [ filefieldprefix int2str(chan_indx(ichan))];
                fieldname = [ filefieldprefix int2str(ichan)];
            else
                fieldname = [ filefieldprefix int2str(subjICs{iSubj}(ichan))];
            end
            all_pac = setfield( all_pac,fieldname , allTrialsPAC{ichan});
        end
        
        % Save  PAC file for each subject
        all_pac.parameters = [];
        for iparam =1:2:length(g.pacparams)
            all_pac.parameters.(g.pacparams{iparam}) = g.pacparams{iparam+1};
        end
        all_pac.parameters    = g.pacparams; % Just the input values. Defaults not printed here
        all_pac.tfparameters  = struct2args(paramstf);
        all_pac.freqs1        = realfreqs1;
        all_pac.freqs2        = realfreqs2;
        all_pac.times         = timesout;
        all_pac.method        = g.method;
        all_pac.trialinfo     = trialinfo;
        all_pac.datatrials    = g.trialindices;
        all_pac.datatype      = datatype;
        all_pac.datafiles     = computeFullFileName( { ALLEEG(inds).filepath }, { ALLEEG(inds).filename });
        if strcmp(prefix,'dat')
             all_pac.labels = dataindx;
        end
        
        if strcmpi(g.savefile, 'on')
            std_savedat( filebase , all_pac );
        end
    end
    
    % Updating STUDY structure
    STUDY.etc.eegpac.method                 = g.method;
    STUDY.etc.eegpac.dataindx               = datindxint;
    STUDY.etc.eegpac.labels                 = dataindx;
    STUDY.etc.eegpac.datatype               = chanorcomp;
    STUDY.etc.eegpac.pactype                = datatype;
    STUDY.etc.eegpac.params.freqs_phase     = all_pac.freqs1;
    STUDY.etc.eegpac.params.freqs_amp       = all_pac.freqs2  ;
    STUDY.etc.eegpac.cache                  = '';
    
    STUDY = pop_comodpacparams(STUDY, 'default');
    STUDY = pop_comodtpacparams(STUDY, 'default');
    STUDY = pop_tfpacparams(STUDY, 'default');
    STUDY = pop_trialspacparams(STUDY, 'default');
    
else
    try
        nsg_info;  % get information on where to create the temporary file
    catch
        error('Plugin nsgportal needs to be in the MATLAB path');
    end
    
    %  Section 1: Create temporary folder and save data
    tmpJobPath = fullfile(outputfolder, 'pactmp');
    if exist(tmpJobPath,'dir'), rmdir(tmpJobPath,'s'); end
    mkdir(tmpJobPath);
    
    % Save data in temporary folder previously created.
    % Here you may change the file name to match the one in the script you will run in NSG
    pop_savestudy( STUDY, ALLEEG, 'filename', EEG.filename, 'filepath', tmpJobPath );
    
    % Copy toolbox to folder. temporary until updated in NSG
    pactoolfolder = which('pop_pacprecomp.m');
    pacpath = fileparts(pactoolfolder);
    copyfile(pacpath,tmpJobPath);
    
    %  Section 2
    %  Manage m-file to be executed in NSG
    %  Write m-file to be run in NSG.
    %  Options defined in plugin are written into the file
    pacnsgnsg_writejobfile
    
    % Section 3
    % Submit job to NSG
    pop_nsg('run',tmpJobPath,'filename', 'pacnsg_job.m', 'jobid', g.jobid,'runtime', g.runtime);
    display([char(10) 'PACTools job (jobID:'  g.jobid ') has been submitted to NSG' char(10) ...
        'Copy or keep in mind the jobID assigned to this job to retreive the results later on.' char(10)...
        'You may follow the status of your job through pop_nsg'...
        char(10)]);
    rmdir(tmpJobPath,'s');
    return;
end

% ------------------------------------|
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

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames)
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end