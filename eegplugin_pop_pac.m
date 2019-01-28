function vers = eegplugin_pop_pac(fig, trystrs, catchstrs);
vers = 'beta02';
plotmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( plotmenu, 'Label', 'ERPAC Tool', 'separator', 'on');
uimenu( submenu, 'label', 'Estim. PAC','callback', 'EEG=pop_pac(EEG); eeglab redraw');
uimenu( submenu, 'label', 'Visualize PAC','callback', 'EEG=pop_plotpac(EEG); eeglab redraw');