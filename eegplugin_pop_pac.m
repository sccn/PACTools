function vers = eegplugin_pop_pac(fig, trystrs, catchstrs);
vers = 'beta03';
plotmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( plotmenu, 'Label', 'ERPAC Tool', 'separator', 'on');
uimenu( submenu, 'label', 'Estimate PAC','callback', 'EEG=pop_pac(EEG); eeglab redraw');
uimenu( submenu, 'label', 'Visualize PAC','callback', 'EEG=pop_plotpac(EEG); eeglab redraw');