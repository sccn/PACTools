function vers = eegplugin_pop_pac(fig,  try_strings, catch_strings);
vers = 'beta03';
plotmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( plotmenu, 'Label', 'ERPAC Tool', 'separator', 'on');
uimenu( submenu, 'label', 'Estimate PAC','callback',...
    [try_strings.no_check '[EEG LASTCOM]=pop_pac(EEG);' catch_strings.add_to_hist]);
uimenu( submenu, 'label', 'Visualize PAC','callback',...
    [try_strings.no_check '[~,LASTCOM]=pop_plotpac(EEG);' catch_strings.add_to_hist]);