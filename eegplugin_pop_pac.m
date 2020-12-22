function vers = eegplugin_pop_pac(fig,  try_strings, catch_strings);
vers = 'beta03';

% Single subject menu
plotmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( plotmenu, 'Label', 'PAC Tools', 'separator', 'on');
uimenu( submenu, 'label', 'Estimate PAC','callback',...
    [try_strings.no_check '[EEG LASTCOM]=pop_pac(EEG);' catch_strings.add_to_hist]);
uimenu( submenu, 'label', 'Plot PAC','callback',...
    [try_strings.no_check '[EEG,LASTCOM]=pop_plotpac(EEG);' catch_strings.add_to_hist]);

% STUDY menu
studymenu = findobj(fig, 'tag', 'study');
submenustudy = uimenu( studymenu, 'Label', 'STUDY PAC Tools', 'separator', 'on','userdata', 'startup:off;study:on');
uimenu( submenustudy, 'label', 'Precompute PAC','callback',[try_strings.no_check '[STUDY, ALLEEG, LASTCOM] = pop_pacprecomp(STUDY, ALLEEG);' catch_strings.add_to_hist],'userdata', 'startup:off;study:on');
uimenu( submenustudy, 'label', 'Plot PAC','callback',[try_strings.no_check '[STUDY, LASTCOM] = pop_stdplotpac(STUDY, ALLEEG);' catch_strings.add_to_hist],'userdata', 'startup:off;study:on');