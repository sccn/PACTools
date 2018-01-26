%% 
% Time specifications:
   Fs = 8000;                   % samples per second
   dt = 1/Fs;                   % seconds per sample
   StopTime = 1;             % seconds
   t = (0:dt:StopTime-dt)';     % seconds
% Sine waves
   figure('Units','normalized','Position',[0.1651 0.1704 0.3573 0.7333]);
   subplot(5,1,1);
   f1 = 9;
   x1 = sin(2*pi*f1*t);
   plot(t,x1,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
   subplot(5,1,2);
   f2 = 5;
   x2 = sin(2*pi*f2*t);
   plot(t,x2,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
   subplot(5,1,3);
   f3 = 7;
   x3 = sin(2*pi*f3*t);
   plot(t,x3,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
   subplot(5,1,4);
   f4 = 15;
   x4 = sin(2*pi*f4*t);
   plot(t,x4,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
    subplot(5,1,5);
   f5 = 60;
   x5 = sin(2*pi*f5*t);
   plot(t,x5,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
   figure;
   X = x1+x5;
   plot(t,X,'Linewidth',3,'Color',[0 0.4470 0.7410]);
   set(gca,'visible','off');
   
   
%%
Amod=1;
Acarr=5;
fmod  = 3;
fcarr = 35;
time = 0:0.001:1;
m=0.5;
sm=Amod*sin(2*pi*fmod*time);
sc=Acarr*sin(2*pi*fcarr*time);
smod=(Amod+m*sm).*sin(2*pi*fcarr*time);
subplot(3,1,1),plot(time,sm,'Linewidth',2,'Color',[0 0.4470 0.7410]);
title('MESSAGE SIGNAL');
set(gca,'visible','off');

subplot(3,1,2),plot(time,sc,'Linewidth',2,'Color',[0 0.4470 0.7410])
title('CARRIER SIGNAL');
set(gca,'visible','off');

subplot(3,1,3)
plot(time,smod,'Linewidth',2,'Color',[0 0.4470 0.7410]);
hold on;
[hAx,hLine1,hLine2] = plotyy(time,abs(hilbert(smod)),time,angle(hilbert(smod)));
hLine1.Color = [1 0 0]; hLine1.LineWidth = 3;
hLine2.Color = [0 1 0]; hLine2.LineWidth = 2; hAx(2).YLim = [-10 15];
title('AMPLITUDE MODULATED SIGNAL');
set(hAx(1),'visible','off');
set(hAx(2),'visible','off');
hlegend = legend('AM Signal','Amplitude','Phase');
hlegend.FontSize = 15;
set(hlegend,'box', 'on','EdgeColor','None');
set(get(hlegend,'BoxFace'),'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.6]));
%%
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','eeglab_data.set','filepath','/Users/amon-ra/WORK/data_sample_eeglab/toplay/sample_data/');

figure('Units','Normalized','Position',[-0.6807 0.5481 0.6906 0.2407]); 
hplot1 = plot(EEG.times,EEG.data(10,:));
axis tight
set(gca,'visible','off');

EEGband1 = pop_eegfiltnew(EEG, 35, 60, 106, 0, [], 1);
figure('Units','Normalized','Position',[-0.6807 0.5481 0.6906 0.2407]); 
hplot2 = plot(EEGband1.times(1000:1500),EEGband1.data(10,1000:1500));
axis tight
set(gca,'visible','off');

figure('Units','Normalized','Position',[-0.6807 0.5481 0.6906 0.2407]); 
% hplot2 = plot(EEGband1.times(1000:1500),EEGband1.data(10,1000:1500),'color',[0.6 0.8 1]); hold on;
axis tight
 plot(EEGband1.times(1000:1500),abs(hilbert(EEGband1.data(10,1000:1500))),'color',[1 0 0],'Linewidth',3); hold on;

set(gca,'visible','off');






