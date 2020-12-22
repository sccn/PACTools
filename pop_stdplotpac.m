% pop_stdpacplot() - Call GUI to compute cross-frequency-coupling coupling.
%             Second level function to compute CFC by calling eeg_pac.m
% Usage:
%   >>  STUDY = pop_stdplotpac(STUDY, ALLEEG);   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG
%   datasets in ALLEEG.
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user edits,
%                if any. Plotted channel measure means (maps, ERSPs, etc.) are added to 
%                the STUDY structure after they are first plotted to allow quick replotting.  
%   com        - Command executed (on development)
%
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2019
%
% Copyright (C) 2019  Ramon Martinez-Cancino,INC, SCCN
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

function [STUDY, com] = pop_stdplotpac(varargin)

com = []; STUDY = [];
icadefs;

if nargin < 2
    help pop_stdplotpac;
    return;
end

%Fixed vals for both entries
AllMethod_listgui = {'Mean vector length modulation index (Canolty et al., 2006)',...
                     'Kullback-Leibler modulation index (Tort et al., 2010)',...
                     'General linear model (Penny et al., 2008)',...
                     'Phase Locking Value (Lachaux et al., 1999)',...
                     'Instantaneous MIPAC (Martinez-Cancino et al., 2019)',...
                     'Event related MIPAC (Martinez-Cancino et al., 2019)'};

AllMethods =  {'mvlmi','klmi','glm','plv','instmipac', 'ermipac'};
alldatatypes = {'Channels', 'ICs Clusters'};
alldatatypes_list = {'channels', 'clusters'};

if ~ischar(varargin{1})  
     STUDY = varargin{1};
     ALLEEG = varargin{2};
     
     % Check for fields here eegpac and pacplotopt
     if length(STUDY)>1
          error('pop_stdplotpac(): Invalid STUDY input');
     end
    if ~isfield(STUDY.etc, 'eegpac') || isempty(STUDY.etc.eegpac) || ~isfield(STUDY.etc, 'pacplotopt')
        error('pop_stdplotpac(): PAC has not been computed for this dataset');
    end
    
    datatype = STUDY.etc.eegpac.datatype; % get data type
         
    % Callbacks
    cb_comod_opt        = 'pop_stdplotpac(''comod_opt'', gcf);';
    cb_comodt_opt       = 'pop_stdplotpac(''comodt_opt'', gcf);';
    cb_tfpac_opt        = 'pop_stdplotpac(''tfpac_opt'', gcf);';
    cb_trialpac_opt     = 'pop_stdplotpac(''trialbasedpac_opt'', gcf);';
    cb_comodplot_clust  = 'pop_stdplotpac(''comodclust'', gcf);';
    cb_comodplot_comp   = 'pop_stdplotpac(''comodcomp'', gcf);';    
    cb_comodtplot_clust = 'pop_stdplotpac(''comodtclust'', gcf);';
    cb_comodtplot_comp  = 'pop_stdplotpac(''comodtcomp'', gcf);';
    cb_tfplot_clust     = 'pop_stdplotpac(''tfpacclust'', gcf);';
    cb_tfplot_comp      = 'pop_stdplotpac(''tfpaccomp'', gcf);';
    cb_trialplot_clust  = 'pop_stdplotpac(''trialpacclust'', gcf);';
    cb_trialplot_comp   = 'pop_stdplotpac(''trialpaccomp'', gcf);';
    cb_show_clust       = 'pop_stdplotpac(''showclust'', gcf);';
    cb_show_comps       = 'pop_stdplotpac(''showcomplist'', gcf);';
    cb_stat_opt         = 'pop_stdplotpac(''stat_opt'',gcf);';
 
    % Gui below    
    guimethodlist     = AllMethod_listgui(find(ismember(AllMethods,STUDY.etc.eegpac.method)));
    
    ph_freqrange  = ['[' num2str(minmax(STUDY.etc.eegpac.params.freqs_phase)) ']'];
    amp_freqrange = ['[' num2str(minmax(STUDY.etc.eegpac.params.freqs_amp)) ']'];
    ph_nfreqs = num2str(numel(STUDY.etc.eegpac.params.freqs_phase));
    amp_nfreqs = num2str(numel(STUDY.etc.eegpac.params.freqs_amp));
    
    comod_enable    = 'on';
    comodt_enable   = 'on';
    freqtime_enable = 'on';
    trialpac_enable = 'on';
    
    % userdata below
    % --------------
    fig_arg{1}{1} = STUDY;
    fig_arg{1}{2} = ALLEEG;
    
    if datatype==1
         labellist1 = 'Select channel to plot';
        labellist2 = 'Select subject to plot';
    else
        labellist1 = 'Select cluster to plot';
        labellist2 = 'Select component to plot        ';
    end
    
    uilist   = { ...
        {'style' 'text' 'string' 'Data type: '  'FontWeight' 'Bold' }      {'style' 'text' 'string' alldatatypes{datatype}} ...
        {'style' 'text' 'string' 'PAC Method: ' 'FontWeight' 'Bold' }      {'style' 'text' 'string' guimethodlist} ...
                                                                           {'style' 'text' 'string' 'Freq range [lo hi] (Hz)' 'FontWeight' 'bold' } {'style' 'text' 'string' '# Frequencies' 'fontweight' 'bold'}...
        {'style' 'text' 'string' 'Phase data'   'FontWeight' 'Bold' }   {'style' 'text' 'string' ph_freqrange}                                   {'style' 'text' 'string'  ph_nfreqs} ...
        {'style' 'text' 'string' 'Amplitude data ' 'FontWeight' 'Bold' }   {'style' 'text' 'string' amp_freqrange}                                  {'style' 'text' 'string'  amp_nfreqs} ...
        {'style' 'text'       'string' 'Select design:' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} {'style' 'popupmenu'  'string' { STUDY.design.name } 'FontWeight' 'Bold' 'tag' 'design' 'value' STUDY.currentdesign } ...
        {'style' 'text'       'string' labellist1 'FontWeight' 'Bold' } {} {'style' 'text'       'string' labellist2 'FontWeight' 'Bold'} ...
        {'style' 'listbox'    'string'  '' 'value' 1 'tag' 'clust_list' 'Callback' cb_show_comps } {'style' 'pushbutton' 'enable'   'on'       'string' 'STATS' 'Callback' cb_stat_opt }     {'style' 'listbox'    'string' '' 'tag' 'comp_list' 'max' 2 'min' 1} ... 
        {'style' 'pushbutton' 'enable'  comod_enable    'string' 'Comodulogram (freqs x freqs)'              'Callback' cb_comodplot_clust   'tag' 'plot_comod_clust'}    {'style' 'pushbutton' 'enable'   comod_enable   'string' 'Params' 'Callback' cb_comod_opt    'tag' 'comod_opt'}     {'style' 'pushbutton' 'enable'  comod_enable    'string' 'Comodulogram (freqs x freqs)'              'Callback' cb_comodplot_comp   'tag' 'plot_comod_comp'} ...
        {'style' 'pushbutton' 'enable'  comodt_enable   'string' 'Time-comodulogram (times x freqs x freqs)' 'Callback' cb_comodtplot_clust 'tag' 'plot_comodt_clust'}   {'style' 'pushbutton' 'enable'   comodt_enable   'string' 'Params' 'Callback' cb_comodt_opt   'tag' 'comodt_opt'}    {'style' 'pushbutton' 'enable'  comodt_enable   'string' 'Time-comodulogram (times x freqs x freqs)' 'Callback' cb_comodtplot_comp 'tag' 'plot_comodt_comp'} ...
        {'style' 'pushbutton' 'enable'  freqtime_enable 'string' 'Time-frequency PAC (freqs x times)'        'Callback' cb_tfplot_clust     'tag' 'plot_tf_clust'}       {'style' 'pushbutton' 'enable'   freqtime_enable 'string' 'Params' 'Callback' cb_tfpac_opt    'tag' 'tfpac_opt'}     {'style' 'pushbutton' 'enable'  freqtime_enable 'string' 'Time-frequency PAC (freqs x times)'        'Callback' cb_tfplot_comp     'tag' 'plot_tf_comp'}     ...
        {'style' 'pushbutton' 'enable'  trialpac_enable 'string' 'Trial-based PAC (Trials x times)'          'Callback' cb_trialplot_clust  'tag' 'plot_trialpac_clust'} {'style' 'pushbutton' 'enable'   trialpac_enable 'string' 'Params' 'Callback' cb_trialpac_opt 'tag' 'trialpac_opt'}  {'style' 'pushbutton' 'enable'  trialpac_enable 'string' 'Trial-based PAC (Trials x times)'          'Callback' cb_trialplot_comp  'tag' 'plot_trialpac_comp'}...
        {}};
    
    guiheight = 14;
    guiwidth = 2.05;
    %%
    geometry    = {{guiwidth guiheight [0 0]    [1 1]}    {guiwidth guiheight [0.4 0]   [1  1]}...
                   {guiwidth guiheight [0 1]    [1 1]}    {guiwidth guiheight [0.4 1]   [1  1]}...
                                                          {guiwidth guiheight [0.6 2]   [1 1]}      {guiwidth guiheight [1.4 2] [1  1]}...
                   {guiwidth guiheight [0 3]   [1 1]}     {guiwidth guiheight [0.62 3]   [1  1]}    {guiwidth guiheight [1.5 3] [1  1]}...
                   {guiwidth guiheight [0 4]   [1 1]}     {guiwidth guiheight [0.62 4]   [1  1]}     {guiwidth guiheight [1.5 4] [1  1]}...                
                   {guiwidth guiheight [0 5]   [1 1]}     {guiwidth guiheight [0.4 5]   [1.65  1]}...
                   {guiwidth guiheight [0 6]   [0.9 1]}   {guiwidth guiheight [0.4 6]   [1 1]}      {guiwidth guiheight [1.15 6]    [0.9 1]}...
                   {guiwidth guiheight [0 7]   [0.9 3]}  {guiwidth guiheight [0.9 7]  [0.25 3]}   {guiwidth guiheight [1.15 7]     [0.9 3]}...
                   {guiwidth guiheight [0 10]  [0.9 1]}  {guiwidth guiheight [0.9 10] [0.25 1]}   {guiwidth guiheight [1.15 10]    [0.9 1]}...
                   {guiwidth guiheight [0 11]  [0.9 1]}  {guiwidth guiheight [0.9 11] [0.25 1]}   {guiwidth guiheight [1.15 11]    [0.9 1]}...
                   {guiwidth guiheight [0 12]  [0.9 1]}  {guiwidth guiheight [0.9 12] [0.25 1]}   {guiwidth guiheight [1.15  12]   [0.9 1]}...
                   {guiwidth guiheight [0 13]  [0.9 1]}  {guiwidth guiheight [0.9 13] [0.25 1]}   {guiwidth guiheight [1.15  13]   [0.9 1]}...
                   {guiwidth guiheight [0 13]  [0.9 1]} };
    
    [out_param userdat tmp res] = inputgui( 'geom' , geometry, 'uilist', uilist,'helpcom', 'pophelp(''pop_stdplotpac'')',...
                                            'title', 'Plot PAC for STUDY -- pop_stdplotpac()' , 'userdata', fig_arg, 'eval', cb_show_clust);
                                        
    if ~isempty(userdat)
        STUDY = userdat{1}{1};
    end
    
else
    hdl = varargin{2};  %figure handle
    userdat  = get(varargin{2}, 'userdat');
    STUDY   = userdat{1}{1};
    ALLEEG  = userdat{1}{2};
    datatype = STUDY.etc.eegpac.datatype;
        
    % Check cluster/channels selected
    ListBoxClustersObj    = findobj(hdl,'tag', 'clust_list');
    ListBoxClustersObjVal = get(ListBoxClustersObj,'value');
    SelectedClustName = STUDY.etc.eegpac.labels{ListBoxClustersObjVal};
    
    ListBoxCompObj    = findobj(hdl,'tag', 'comp_list');
    ListBoxCompObjVal = get(ListBoxCompObj,'value');
    
    % Components or subjects
    if datatype==1
        if ListBoxCompObjVal(1) ~= 1
            %         channels
            subject = STUDY.design(STUDY.currentdesign).cases.value{ListBoxCompObjVal-1};
        end
    else
        clsindx = find(~cellfun('isempty',strfind({STUDY.cluster.name},STUDY.etc.eegpac.labels{ListBoxClustersObjVal})));
        subcomplist = {STUDY.datasetinfo(STUDY.cluster(clsindx).sets(1,:)).subject};
        allcomplist = STUDY.cluster(clsindx).comps;
    end
    
%     if datatype == 1
%         SelectedCompName  = STUDY.etc.eegpac.labels{ListBoxClustersObjVal};
%     else
%        SelectedCompName = [];
%     end
    SelectedCompName  = STUDY.etc.eegpac.labels{ListBoxClustersObjVal};
    design  = get(findobj('parent', hdl, 'tag', 'design')      , 'value');
            
    try
        switch  varargin{1}
            case {'comodclust', 'comodtclust', 'tfpacclust','trialpacclust'}
                
                plotting_option = varargin{1};
                plotting_option = [ 'plot' plotting_option(1:end-5)];
                if datatype == 1
                    a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype} ''',{'  vararg2str(SelectedClustName) '}, ''design'', ' int2str(design) ');' ];
                else
                    a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype} ''','  num2str(clsindx) ', ''design'', ' int2str(design) ');' ];
                    
                end
                eval(a); 
                %STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                
                 % Update Study history
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat);

           case {'comodcomp', 'comodtcomp', 'tfpaccomp','trialpaccomp'}
               
                plotting_option = varargin{1};
                plotting_option = [ 'plot' plotting_option(1:end-4)];     
                
                if ListBoxCompObjVal(1) ~= 1  % check that not all onechan in channel are requested
                    if datatype == 1
                     a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype}  ''',{'  vararg2str({SelectedCompName}) '}, ''subject'', ''' subject ''', ''design'', ' int2str(design) ' );' ];
                    else   
                     a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype}  ''','  num2str(clsindx) ',''comps'','  num2str(ListBoxCompObjVal-1)  ', ''design'', ' int2str(design) ' );' ];
                    end
                     eval(a); 
                     %STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                else
                    if datatype == 1
                        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype}  ''',{'  vararg2str({SelectedCompName}) '}, ''plotsubjects'', ''on'', ''design'', ' int2str(design) ' );' ];
                    else
                        a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''' alldatatypes_list{datatype}  ''','  num2str(clsindx) ', ''plotsubjects'', ''on'', ''design'', ' int2str(design) ' );' ];
                    end
                    eval(a); 
                    %STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                end
                 
                % Update Study history
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat); 
                
            case  'comod_opt'
                STUDY = pop_comodpacparams(STUDY);
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat);
                
            case  'comodt_opt'
                STUDY = pop_comodtpacparams(STUDY);
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat);
                
            case  'tfpac_opt'
                STUDY = pop_tfpacparams(STUDY);
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat);
                
            case  'trialbasedpac_opt'
                STUDY = pop_trialspacparams(STUDY);
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat);
                
            case  'enablebuttons'      
                 switch STUDY.etc.eegpac.pactype
                     case 'PAC'
                         set(findobj(hdl,'tag', 'plot_comod_clust'), 'enable', 'on');      set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'on');      set(findobj(hdl,'tag', 'plot_comod_comp'), 'enable', 'on');  
                         set(findobj(hdl,'tag', 'plot_comodt_clust'), 'enable', 'on');     set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'on');     set(findobj(hdl,'tag', 'plot_comodt_comp'), 'enable', 'on');  
                         set(findobj(hdl,'tag', 'plot_tf_clust'), 'enable', 'on');         set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'on');      set(findobj(hdl,'tag', 'plot_tf_clust_comp'), 'enable', 'on');    
                         set(findobj(hdl,'tag', 'plot_trialpac_clust'), 'enable', 'off');  set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'off');  set(findobj(hdl,'tag', 'plot_trialpac_comp'), 'enable', 'off');
                     case 'MIPAC'
                         set(findobj(hdl,'tag', 'plot_comod_clust'), 'enable', 'on');      set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'on');      set(findobj(hdl,'tag', 'plot_comod_comp'), 'enable', 'on');  
                         set(findobj(hdl,'tag', 'plot_comodt_clust'), 'enable', 'on');     set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'on');     set(findobj(hdl,'tag', 'plot_comodt_comp'), 'enable', 'on');  
                         set(findobj(hdl,'tag', 'plot_tf_clust'), 'enable', 'on');         set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'on');      set(findobj(hdl,'tag', 'plot_tf_clust_comp'), 'enable', 'on');    
                         set(findobj(hdl,'tag', 'plot_trialpac_clust'), 'enable', 'on');  set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'off');  set(findobj(hdl,'tag', 'plot_trialpac_comp'), 'enable', 'on');
                 end    
            case 'showclust'
                if datatype == 1
                    for i =1: length(STUDY.etc.eegpac.labels)
                        displayclust{i} = ['All ' STUDY.etc.eegpac.labels{i}];  
                    end
                    set(findobj(hdl,'tag', 'clust_list'),'String',displayclust, 'Value', 1);
                else
                    for i =1: length(STUDY.etc.eegpac.labels)
                        clsindx = find(~cellfun('isempty',strfind({STUDY.cluster.name},STUDY.etc.eegpac.labels{i})))
                        displayclust{i} = [STUDY.etc.eegpac.labels{i} ' (' num2str(length(STUDY.cluster(clsindx).comps)) ' ICs)' ];  
                    end
                    set(findobj(hdl,'tag', 'clust_list'),'String',displayclust, 'Value', 1);
                end
                
                pop_stdplotpac('showcomplist', hdl);  % Show componenets/subjects in list
                pop_stdplotpac('enablebuttons', hdl); % Enable/disable buttons based on metric
                
            case 'showcomplist'
                 if datatype == 1
                     allsubjects = STUDY.design(STUDY.currentdesign).cases.value;
                     complist = {'All subjects'};
                    for l = 1:length(allsubjects)
                        complist{end+1} = [ allsubjects{l} ' ' SelectedClustName ];
                    end
                else
                     clsindx = find(~cellfun('isempty',strfind({STUDY.cluster.name},STUDY.etc.eegpac.labels{ListBoxClustersObjVal})));
                     subcomplist = {STUDY.datasetinfo(STUDY.cluster(clsindx).sets(1,:)).subject};
                     allcomplist = STUDY.cluster(clsindx).comps;
                     
                     complist = {'All components'};
                    for l = 1:length(allcomplist)
                        complist{end+1} = [ subcomplist{l} ' IC' num2str(allcomplist(l)) ];
                    end
                 end
                set(findobj(hdl,'tag', 'comp_list'),'String',complist, 'Value', 1);
                
            case 'stat_opt' % save the list of selected channels
                [STUDY com] = pop_statparams(STUDY);
                userdat{1}{1} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)
        end
    catch
        eeglab_error;
    end
end