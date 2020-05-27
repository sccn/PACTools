% pop_genpac() - Generate an EEG set with a syntethic PAC signal. 
%                The modulation is inserted in alternated segments. 
%                The waveform  of the PAC segments can be defined as well.
%                         
% Usage:
%   >> [EEG, data_pac, time] = pop_genpac(fc,fm,tlimits);
%   >> [EEG, data_pac, time] = pop_genpac(fc,fm,tlimits,'srate',srate,'cpfunc','block','blockamp',1,'plot_flag',1,'padtime',padtime,'snr',snrval);
% Inputs:
%    fc       - Frequency of the carrier [Hz]
%    fm       - Frequency of the modulator. [Hz]
%    tlimits  - Minimun and maximun time for simulation in seconds (s). [min max]
%
% Optional inputs
%   cpfunc    - {'linear','sine', 'exp', 'block'} Function to use to generate the waveform in the segments
%               with PAC ~= 0. Default {'linear'}
%   srate     - Sampling frequency [Hz]. Default [500]
%   Ac        - Maximun amplitude of the carrier signal. Default [5] 
%   Am        - Maximun amplitude of the modulator signal . Default [1] 
%   snr       - Signal to noise ratio for the simulated signal. Default Inf
%   nsegm     - Number of segments to breake the signal into. This is to 
%               create alternated segments with PAC and no PAC. default [3]
%  linslope   - Slope of the line if cpfunc = {'linear'} Default [1]
%  sinamp     - Sine max amplitude if cpfunc = 'sine'. Default [1]
%  expamp     - Exponential max amplitude if cpfunc = 'exp'. Default [1]
%  blockamp   - Block max amplitude if cpfunc = 'block'. Default [1]
%  fsin       - Frequency of the sin function if cpfunc = 'sine'. Default [0.1]
%  padtime    - Time padding with 0 at the beggining and the end of the simulated
%               signal. Default [0]
% plot_flag   - Plot signal generated [0,1]. Default : no plot[0]
% 
% Outputs: 
% amplitude_mod     - Amplitude modulated signal
% phase_signal      - Time serie of the phase of the cmodulator used to generate the signal
% time                 - Time points of the generated signal
% m                 - Time serie of the modulation strength
%
%
% Author: Ramon Martinez-Cancino, SCCN/INC, UCSD 2018
%
%
% Author: Ramon Martinez-Cancino, SCCN, 2018
%
% Copyright (C) 2018  Ramon Martinez-Cancino, 
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
function [EEG, data_pac, time] = pop_genpac(fc,fm,tlimits,varargin)

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('pop_genpac() error: calling convention {''key'', value, ... } error'); return;
end;
defaultname = ['Simpac_famp_' num2str(fc) '_fphs_' num2str(fm) '.set'];
setname     = ['PAC signal famp_' num2str(fc) '_fphs_' num2str(fm)];

try g.plot_flag;   catch, g.plot_flag   = 0;             end;
try g.cpfunc;      catch, g.cpfunc      = 'linear';      end;
try g.srate;       catch, g.srate       = 500;           end;
try g.Ac;          catch, g.Ac          = 5;             end;
try g.Am;          catch, g.Am          = 1;             end;
try g.snr;         catch, g.snr         = Inf;           end;
try g.nsegm;       catch, g.nsegm       = 3;             end;
try g.linslope;    catch, g.linslope    = 1;             end;
try g.sinamp;      catch, g.sinamp      = 1;             end;
try g.expamp;      catch, g.expamp      = 1;             end;
try g.blockamp;    catch, g.blockamp    = 1;             end;
try g.m;           catch, g.m           = 0.5;           end;
try g.fsin;        catch, g.fsin        = 0.1;           end;
try g.padtime;     catch, g.padtime     = 0;             end;
try g.maxshift;    catch, g.maxshift    = 0;             end;
try g.nchan;       catch, g.nchan       = 1;             end;
try g.ntrials;     catch, g.ntrials     = 1;             end;
try g.filename;    catch, g.filename    = defaultname;   end;
try g.setname;     catch, g.setname     = setname;       end;

% GUI here
if nargin < 1
    disp('GUI under development');
end

m = g.m;
if (any(g.m < 0)||any(g.m > 1))
    error('Coupling value should be less than or equal to one and greater than zero');
end
% General settings
time  = tlimits(1):1/g.srate:tlimits(2)+2*g.padtime;  % Total time for simulation
amplitude_mod = zeros(1,length(time));
phase_signal  = zeros(1,length(time));

% Carrier signal generation
carrier_signal = g.Ac*sin(2*pi*fc*time);

%% Coupling functions. Determining shape of m
if ~isempty(g.cpfunc)
    m        = zeros(1,length(time));
    npts_seg = floor((length(time)-2*g.padtime*g.srate)/g.nsegm);
    count    = 1;
    treset   = 0;
    if g.padtime~=0
        count  = g.padtime*g.srate + 1;
        treset = g.padtime;
    end
    for iseg = 1:g.nsegm
        if iseg ~= g.nsegm
            ti = time(count:count + npts_seg-1);
            switch lower(g.cpfunc)
                case {'linear'}
                    m(count:count + npts_seg-1) = g.linslope*(ti-treset);
                case {'linearsegm'}
                    if rem(iseg,2) ~= 0, mi = 0; else mi = g.linslope; end
                    m(count:count + npts_seg-1) = mi*(ti-treset);
                case 'sine'
                    m(count:count + npts_seg-1) = abs(g.sinamp*sin(2*pi*g.fsin*(ti-treset)));    
                case 'exp'
                    if rem(iseg,2) ~= 0, mi = 0; else mi = g.expamp; end
                    m(count:count + npts_seg-1) = mi*(exp(ti-treset)-1);    
                case 'block'
                    if rem(iseg,2) ~= 0, mi = 0; else mi = g.blockamp; end
                    m(count:count + npts_seg-1) = mi;
            end
            count = count + npts_seg;
            treset = ti(end);
        else
            if g.padtime==0
                ti = time(count:end);
            else
                ti = time(count:end-g.padtime*g.srate);
            end
            switch lower(g.cpfunc)
                case {'linear'}
                    m(count:count + length(ti)-1) = g.linslope*(ti-treset);
                    case 'linearsegm'
                    if mi == 0, mi = g.linslope; else mi = 0; end
                    m(count:count + length(ti)-1) = mi*(ti-treset);
                     m = m+abs(min(m));
                case 'sine'
                    m(count:count + length(ti)-1) = abs(g.sinamp*sin(2*pi*g.fsin*(ti-treset)));
                case 'exp'
                    if mi == 0, mi = g.blockamp; else mi = 0; end
                    m(count:count + length(ti)-1) = mi*(exp(ti-treset)-1);
                    m = m+abs(min(m));
                case 'block'
                    if mi == 0, mi = g.blockamp; else mi = 0; end
                    m(count:count + length(ti)-1) = mi;
            end
        end
    end  
end
mmod = m./max(m);

% Amplitude Modulation
if length(mmod) == length(amplitude_mod)
    ti = time;
    amplitude_mod  = (1+ mmod.*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod.*cos(2*pi*fm*ti);
    phase_signal   = g.Am*cos(2*pi*fm*ti);
    amplitude_mod  = awgn(amplitude_mod,g.snr,'measured');
else
    npts_seg = floor(length(time)/length(mmod));
    count = 1;
    
    for iseg = 1: length(mmod)
        if iseg ~= length(mmod)
            ti = time(count:count + npts_seg-1);
            amplitude_mod(count:count + npts_seg-1)  = (1+ mmod(iseg).*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod(iseg).*cos(2*pi*fm*ti);
            amplitude_mod(count:count + npts_seg-1) = awgn(amplitude_mod(count:count + npts_seg-1),g.snr,'measured');
            % Modulating signal generation
            phase_signal(count:count + npts_seg-1) = cos(2*pi*fm*ti);
            count = count + npts_seg;
        else
            ti = time(count:end);
            phase_signal(count:count+ length(ti)-1) = cos(2*pi*fm*ti);
            amplitude_mod(count:count+length(ti)-1) = (1+ mmod(iseg).*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod(iseg).*cos(2*pi*fm*ti);
            amplitude_mod(count:count+length(ti)-1) = awgn(amplitude_mod(count:count+length(ti)-1), g.snr,'measured');
        end
    end
end

% Plots
if g.plot_flag
    figure('Units','normalized','Position',[0.3049 0.2544 0.5486 0.6311]);
    ax(1) = subplot(3,1,1);
    plot(time,phase_signal);
    title (['Modulator/Phase f = ' num2str(fm) 'Hz']);
    ylabel ('Amplitude');
    ax(1).XTickLabel = [];
    grid on;
    
    % Plot 2. Carrier
    ax(2) = subplot(3,1,2);
    plot(time,carrier_signal);
    title (['Carrier/Amplitude f = ' num2str(fc) 'Hz']);
    ylabel ('Amplitude');
    ax(2).XTickLabel = [];
    grid on;
    
    % Plot 3. Modulated signal
    ax(3) =subplot(3,1,3);
    plot(time,amplitude_mod); hold on;
    plot(time,m*max(amplitude_mod),'color',[215 25 25]./255,'linewidth',2)
    title ('Amplitude Modulated Signal');
    xlabel ('time(sec)');
    ylabel ('Amplitude');
    grid on;
    hlegend = legend('Signal', 'Modulation');
    hlegend.Box = 'off';
    hlegend. Position = [0.755 0.2650 0.1032 0.0431];
    linkaxes(ax,'x');
    
    % Format data
    if g.maxshift~=0
        s = randi(g.maxshift, g.ntrials,1);
    end
   for ichan = 1:g.nchan
    for itrial=1:g.ntrials
        if g.snr ~= Inf
        datatmp = awgn(amplitude_mod,g.snr,'measured');
        else
           datatmp = amplitude_mod; 
        end
        if g.maxshift~=0
            datatmp = [datatmp(s(itrial):end) datatmp(1:s(itrial)-1)];
            data_pac(ichan,:,itrial) = [datatmp(end-ceil(g.maxshift/2):end) datatmp(1:end-ceil(g.maxshift/2)-1)];
        else
            data_pac(ichan,:,itrial) = datatmp;
        end  
    end
   end
    
    % Creating EEG structure
    EEG = eeg_emptyset;
    EEG.setname  = g.setname;
    EEG.filename = g.filename;
    EEG.nbchan   = g.nchan;
    EEG.trials   = g.ntrials;
    EEG.pnts     = size(data_pac,2);
    EEG.xmin     = time(1);
    EEG.xmax     = time(end);
    EEG.data     = data_pac;
    EEG.srate    = g.srate;
    EEG.times    = time;
end