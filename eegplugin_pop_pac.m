function vers = eegplugin_pop_pac(fig, trystrs, catchstrs);
vers = 'beta01';
plotmenu = findobj(fig, 'tag', 'tools');
uimenu( plotmenu, 'label', 'pop_pac','callback', 'EEG=pop_pac(EEG); eeglab redraw');