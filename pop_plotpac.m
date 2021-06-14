% pop_plotpac() - Call GUI to compute cross-frequency-coupling coupling.
%             Second level function to compute CFC by calling eeg_pac.m
% Usage:
%   >>  pac = pop_plotpac(EEG);
%
% Inputs:
%  EEG          - [Structure] Input dataset as an EEGLAB EEG structure

% Outputs:
%  EEG          - [Structure] EEG dataset structure with pac results
%  com          - Command executed (on development)
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

function [EEG, com] = pop_plotpac(varargin)
com = []; EEG = [];

%Fixed vals for both entries
AllMethod_listgui = {'Mean vector length modulation index (Canolty et al., 2006)',...
                     'Kullback-Leibler modulation index (Tort et al., 2010)',...
                     'General linear model (Penny et al., 2008)',...
                     'Phase Locking Value (Lachaux et al., 1999)',...
                     'Instantaneous MIPAC (Martinez-Cancino et al., 2019)',...
                     'Event related MIPAC (Martinez-Cancino et al., 2019)'};

AllMethods =  {'mvlmi','klmi','glm','plv','instmipac', 'ermipac'};
alldatatypes = {'Channels', 'Components (ICs)'};

if ~ischar(varargin{1})
     EEG = varargin{1};
     
     % Check for fields here eegpac and pacplotopt
     if length(EEG)>1
          error('pop_plotpac(): Invalid EEG input');
     end
     
    if ~isfield(EEG.etc, 'eegpac') || isempty(EEG.etc.eegpac) || ~isfield(EEG.etc, 'pacplotopt')
        error('pop_plotpac(): PAC has not been computed for this dataset');
    end
     
    % Callbacks
    cb_comod_opt    = ['pop_plotpac(''comod_opt'', gcf);'];
    cb_comodt_opt   = ['pop_plotpac(''comodt_opt'', gcf);'];
    cb_tfpac_opt    = ['pop_plotpac(''tfpac_opt'', gcf);'];
    cb_trialpac_opt = ['pop_plotpac(''trialbasedpac_opt'', gcf);'];
    cb_lboxchanpair = ['pop_plotpac(''enablebuttons'', gcf);'];
    
    cb_comodplot  = ['pop_plotpac(''plot_comod'', gcf);'];
    cb_comodtplot = ['pop_plotpac(''plot_comodt'', gcf);'];
    cb_tfplot     = ['pop_plotpac(''plot_tfpac'', gcf);'];
    cb_trialplot  = ['pop_plotpac(''plot_trialpac'', gcf);'];

    
    % Gui below
    datatype = EEG.etc.eegpac(1).datatype;
    pair_list = [EEG.etc.eegpac.labels]; % If results in cell this mus be changed. hera assuming using structure
    defaultmethodlist = fieldnames(EEG.etc.eegpac(1));
    guimethodlist = AllMethod_listgui(find(ismember(AllMethods,fieldnames(EEG.etc.eegpac(1)))));
    
    ph_freqrange  = ['[' num2str(myminmax(EEG.etc.eegpac(1).params.freqs_phase)) ']'];
    amp_freqrange = ['[' num2str(myminmax(EEG.etc.eegpac(1).params.freqs_amp)) ']'];
    ph_nfreqs = num2str(numel(EEG.etc.eegpac(1).params.freqs_phase));
    amp_nfreqs = num2str(numel(EEG.etc.eegpac(1).params.freqs_amp));
    
    comod_enable    = 'on';
    comodt_enable   = 'on';
    freqtime_enable = 'on';
    trialpac_enable = 'on';
    
    % userdata below
    % --------------
    fig_arg{1}{1} = EEG;
       
    uilist   = { ...
        {'style' 'text' 'string' 'Data type: ' 'FontWeight' 'Bold' }  {'style' 'text' 'string' alldatatypes{datatype}} ...
                                                                           {'style' 'text' 'string' 'Freq range [lo hi] (Hz)' 'FontWeight' 'bold' } {'style' 'text' 'string' '# Frequencies' 'fontweight' 'bold'}...
        {'style' 'text' 'string' 'Phase data' 'FontWeight' 'Bold' }      {'style' 'text' 'string' ph_freqrange}                                   {'style' 'text' 'string'  ph_nfreqs} ...
        {'style' 'text' 'string' 'Amplitude data ' 'FontWeight' 'Bold' }  {'style' 'text' 'string' amp_freqrange}                                  {'style' 'text' 'string'  amp_nfreqs} ...
        {'style' 'text' 'string' 'Select IC|channel pairing' 'FontWeight' 'Bold' }...
        {'style' 'listbox'  'string' pair_list 'tag' 'lbox_pairs' 'max' 1 'min' 1 'tag' 'lbox_chanpair'} ...
        {'style' 'text'       'string' 'Select method to plot ' 'FontWeight' 'Bold'  }...
        {'style' 'popupmenu'  'string' guimethodlist 'tag' 'pupm_method' 'value' 1 'callback' cb_lboxchanpair} ...
        {'style' 'checkbox' 'tag' 'chckbx_signif' 'value' 0 'string' 'Plot significance if computed'}...
        { }...
        {'style' 'text'       'string' 'Plot PAC ' 'FontWeight' 'Bold'  }...
        {'style' 'pushbutton' 'enable'  comod_enable    'string' 'Comodulogram (freqs x freqs)'          'Callback' cb_comodplot  'tag' 'plot_comod'}    {'style' 'pushbutton' 'enable'   comod_enable    'string' 'Params' 'Callback' cb_comod_opt   'tag' 'comod_opt'} ...
        {'style' 'pushbutton' 'enable'  comodt_enable   'string' 'Time-comodulogram (times x freqs x freqs)' 'Callback' cb_comodtplot 'tag' 'plot_comodt'}   {'style' 'pushbutton' 'enable'   comodt_enable   'string' 'Params' 'Callback' cb_comodt_opt   'tag' 'comodt_opt'} ...
        {'style' 'pushbutton' 'enable'  freqtime_enable 'string' 'Time-frequency PAC (freqs x times)'        'Callback' cb_tfplot     'tag' 'plot_tf'}       {'style' 'pushbutton' 'enable'   freqtime_enable 'string' 'Params' 'Callback' cb_tfpac_opt    'tag' 'tfpac_opt'} ...
        {'style' 'pushbutton' 'enable'  trialpac_enable 'string' 'Trial-based PAC (Trials x times)'       'Callback' cb_trialplot  'tag' 'plot_trialpac'} {'style' 'pushbutton' 'enable'   trialpac_enable 'string' 'Params' 'Callback' cb_trialpac_opt 'tag' 'trialpac_opt'} ...
        {}};
    
    guiheight = 17;
    guiwidth = 2.2;
    
    geometry    = {{guiwidth guiheight [0 0]    [1 1]}  {guiwidth guiheight [0.4 0] [1  1]}...
                                                        {guiwidth guiheight [0.6 1] [1 1]}  {guiwidth guiheight [1.5 1] [1  1]}...
        {guiwidth guiheight [0 2] [1 1]}  {guiwidth guiheight [0.8 2] [1  1]} {guiwidth guiheight [1.7 2] [1  1]}...
        {guiwidth guiheight [0 3] [1 1]}  {guiwidth guiheight [0.8 3] [1  1]} {guiwidth guiheight [1.7 3] [1  1]}...
        {guiwidth guiheight [0 4]    [1 1]}...
        {guiwidth guiheight [0 5]    [2.2 3]}...
        {guiwidth guiheight [0 8]    [1 1]}...
        {guiwidth guiheight [0 9]    [2.1 1]}...
        {guiwidth guiheight [0 10]   [1.6 1]}...
        {guiwidth guiheight [0 11]   [1.6 1]}...
        {guiwidth guiheight [0 12]   [1.6 1]}...
        {guiwidth guiheight [0 13]   [1.75 1]}  {guiwidth guiheight [1.75 13]  [0.45 1]}...
        {guiwidth guiheight [0 14]   [1.75 1]}  {guiwidth guiheight [1.75 14]  [0.45 1]}...
        {guiwidth guiheight [0 15]   [1.75 1]}  {guiwidth guiheight [1.75 15]  [0.45 1]}...
        {guiwidth guiheight [0 16]   [1.75 1]}  {guiwidth guiheight [1.75 16]  [0.45 1]}...
        {guiwidth guiheight [0 17]   [1.5 1]}...
        };
    
    [out_param userdat tmp res] = inputgui( 'geom' , geometry, 'uilist', uilist,'helpcom', 'pophelp(''pop_pacplot'')',...
                                            'title', 'Plot PAC for single subject -- pop_pacplot()' , 'userdata', fig_arg,...
                                            'eval', cb_lboxchanpair );
                                        
    if ~isempty(userdat)
        EEG = userdat{1}{1};
    end
    
else
    hdl = varargin{2};  %figure handle
    userdat  = get(varargin{2}, 'userdat');
    EEG   = userdat{1}{1};
    
    flag_stat = get(findobj(hdl,'tag', 'chckbx_signif'),'Value');
    
    % Check pair selected
    ListBoxObj = findobj(hdl,'tag', 'lbox_chanpair');
    PairIndxVal = get(ListBoxObj,'value');
    
    % Check method in structure
    PopupMenuMethod = findobj(hdl,'tag', 'pupm_method');
    CurrMethodString = get(PopupMenuMethod,'string');
    CurrMethodVal = get(PopupMenuMethod,'value');
    CurrMethod  = CurrMethodString{CurrMethodVal};
    CurrMethodCode =  AllMethods{find(~cellfun(@isempty, strfind(AllMethod_listgui,CurrMethod)))};
    
    try
        switch  varargin{1}
            case 'plot_comod'
                eeg_plotcomod(EEG,'plotindx', PairIndxVal,'pacmethod', CurrMethodCode, 'plotsignif', flag_stat);
                
            case 'plot_comodt'
                eeg_plotcomodt(EEG,'plotindx', PairIndxVal,'pacmethod', CurrMethodCode,'plotsignif', flag_stat);
                
            case 'plot_tfpac'
                 eeg_plottfpac(EEG,'plotindx', PairIndxVal,'pacmethod', CurrMethodCode);
                
            case 'plot_trialpac'
                eeg_plotrialpac(EEG,'plotindx', PairIndxVal,'pacmethod', CurrMethodCode);
                
            case  'comod_opt'
                EEG = pop_comodpacparams(EEG);
                userdat{1}{1} = EEG;
                set(hdl, 'userdat',userdat);
                
            case  'comodt_opt'
                EEG = pop_comodtpacparams(EEG);
                userdat{1}{1} = EEG;
                set(hdl, 'userdat',userdat);
                
            case  'tfpac_opt'
                EEG = pop_tfpacparams(EEG);
                userdat{1}{1} = EEG;
                set(hdl, 'userdat',userdat);
                
            case  'trialbasedpac_opt'
                EEG = pop_trialspacparams(EEG);
                userdat{1}{1} = EEG;
                set(hdl, 'userdat',userdat);
                
            case  'enablebuttons'
                % Dimension of computed PAC
                if ~isempty(EEG.etc.eegpac(PairIndxVal).(CurrMethodCode).pacval)
                    pacdim = EEG.etc.eegpac(PairIndxVal).(CurrMethodCode).dim;
                else
                    pacdim = 0;
                end
                
                 switch pacdim
                     case 0
                         set(findobj(hdl,'tag', 'plot_comod'), 'enable', 'off');    set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'off'); 
                         set(findobj(hdl,'tag', 'plot_comodt'), 'enable', 'off');   set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'off'); 
                         set(findobj(hdl,'tag', 'plot_tf'), 'enable', 'off');       set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'off'); 
                         set(findobj(hdl,'tag', 'plot_trialpac'), 'enable', 'off'); set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'off');  
                     case 1
                         set(findobj(hdl,'tag', 'plot_comod'), 'enable', 'on');     set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_comodt'), 'enable', 'off');   set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'off'); 
                         set(findobj(hdl,'tag', 'plot_tf'), 'enable', 'off');       set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'off'); 
                         set(findobj(hdl,'tag', 'plot_trialpac'), 'enable', 'off'); set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'off');  
                     case 2
                         set(findobj(hdl,'tag', 'plot_comod'), 'enable', 'on');      set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_comodt'), 'enable', 'on');     set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_tf'), 'enable', 'on');         set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_trialpac'), 'enable', 'off');  set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'off');  
                     case 3
                          set(findobj(hdl,'tag', 'plot_comod'), 'enable', 'on');   set(findobj(hdl,'tag', 'comod_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_comodt'), 'enable', 'on');   set(findobj(hdl,'tag', 'comodt_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_tf'), 'enable', 'on');       set(findobj(hdl,'tag', 'tfpac_opt'), 'enable', 'on'); 
                         set(findobj(hdl,'tag', 'plot_trialpac'), 'enable', 'on'); set(findobj(hdl,'tag', 'trialpac_opt'), 'enable', 'on');  
                 end    
        end
    catch
        eeglab_error;
    end
end