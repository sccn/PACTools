% generate_pac_signal() - Generate a signal with coupled phase and amplitude (PAC)
%                         in the frequencies defined in the inputs. The PAC
%                         is inserted in alternated segments The waveform 
%                         of the PAC segments can be defined as well.
%                         
% Usage:
%   >> [data,simspk.t_out,simspk.phase_signal,simspk.m]  = generate_pac_signal(fc,fm,tlimits);
%   >> [data, t_out,phase_signal,m]  = generate_pac_signal(fc,fm,tlimits,'fs',srate,'cpfunc','block','blockamp',1,'plot_flag',1,'padtime',padtime,'snr',snrval);
% Inputs:
%    fc       - Frequency of the carrier [Hz]
%    fm       - Frequency of the modulator. [Hz]
%    tlimits  - Minimun and maximun time for simulation in seconds (s). [min max]
%
% Optional inputs
%   cpfunc    - {'linear','sine', 'exp', 'block'} Function to use to generate the waveform in the segments
%               with PAC ~= 0. Default {'linear'}
%   fs        - Sampling frequency [Hz]. Default [500]
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
% t                 - Time points of the generated signal
% m                 - Time serie of the modulation strength
%
%
% Author: Ramon Martinez-Cancino, SCCN/INC, UCSD 2016
%
%
% Author: Ramon Martinez-Cancino, SCCN, 2016
%
% Copyright (C) 2014  Ramon Martinez-Cancino, 
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
function [amplitude_mod, t,phase_signal,m] = generate_pac_signal(fc,fm,tlimits,varargin)


try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('std_clustinfo() error: calling convention {''key'', value, ... } error'); return;
end;

pop_flag = 0;
try g.plot_flag;   catch, g.plot_flag   = 0;             end;
try g.cpfunc;      catch, g.cpfunc      = 'linear';      end;
try g.fs;          catch, g.fs          = 500;           end;
try g.Ac;          catch, g.Ac          = 5;             end;
try g.Am;          catch, g.Am          = 1;             end;
try g.snr;          catch, g.snr        = Inf;           end;

try g.nsegm;       catch, g.nsegm       = 3;             end;
try g.linslope;    catch, g.linslope    = 1;             end;
try g.sinamp;      catch, g.sinamp      = 1;             end;
try g.expamp;      catch, g.expamp      = 1;             end;
try g.blockamp;    catch, g.blockamp    = 1;             end;
try g.m;           catch, g.m           = 0.5;            end;
try g.fsin;        catch, g.fsin        = 0.1;           end;
try g.padtime;     catch, g.padtime     = 0;             end;
try g.method;      catch, g.method      = 'wiki';             end;

if ~isempty(g.m)
    m = g.m;
elseif isempty(g.m) && isempty(g.cpfunc)
    error('Coupling must be defined');
end

if (any(0>g.m)||any(g.m>1))
    error('m may be less than or equal to one and greater than zero');
end
% General settings
t  = tlimits(1):1/g.fs:tlimits(2)+2*g.padtime;  % Total time for simulation
noise_amp = 0;
amplitude_mod = zeros(1,length(t));
phase_signal  = zeros(1,length(t));

% Carrier signal generation
g.Ac = 5; % A Amplitude
g.Am = 1;
carrier_signal = g.Ac*sin(2*pi*fc*t);

%% Coupling functions. Determining shape of m
if ~isempty(g.cpfunc)
    m        = zeros(1,length(t));
    npts_seg = floor((length(t)-2*g.padtime*g.fs)/g.nsegm);
    count    = 1;
    treset   = 0;
    if g.padtime~=0
        count  = g.padtime*g.fs + 1;
        treset = g.padtime;
    end
    for iseg = 1:g.nsegm
        if iseg ~= g.nsegm
            ti = t(count:count + npts_seg-1);
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
                ti = t(count:end);
            else
                ti = t(count:end-g.padtime*g.fs);
            end
            switch lower(g.cpfunc)
                % LINEAR
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

m = m./max(m);
if strcmp(g.method,'tort')
    mmod = -m+1;
else
    mmod = m;
end
% AM Modulation
if length(mmod) ~= 1
    if length(mmod) == length(amplitude_mod)
        ti = t;
        if strcmp(g.method,'tort')
            amplitude_mod  = (g.Ac*((1-mmod).*sin(2*pi*fm*ti)+1+mmod)/2).*sin(2*pi*fc*ti) + g.Am*sin(2*pi*fm*ti) + noise_amp.*rand(1,length(ti));
            phase_signal   = cos(2*pi*fm*ti);
        elseif strcmp(g.method,'wiki')
%             amplitude_mod  = sin(2*pi*fc*ti) + g.Ac.*mmod.*(sin(2*pi*(fm+fc)*ti) - sin(2*pi*(fm-fc)*ti))/2 + mmod.*cos(2*pi*fm*ti);
            
            amplitude_mod  = (1+ mmod.*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod.*cos(2*pi*fm*ti) + noise_amp.*rand(1,length(ti));
            phase_signal   = g.Am*cos(2*pi*fm*ti);
            
            
%             amplitude_mod  = sin(2*pi*fc*ti) + mmod.*sin(2*pi*fm*ti).*sin(2*pi*fc*ti) + 0.1.*sin(2*pi*fm*ti);
%             phase_signal   = sin(2*pi*fm*ti);
        end
        amplitude_mod  = awgn(amplitude_mod,g.snr,'measured');
%         phase_signal   = cos(2*pi*fm*ti);
    else
        npts_seg = floor(length(t)/length(mmod));
        count = 1;
        
        for iseg = 1: length(mmod)
            if iseg ~= length(mmod)
                ti = t(count:count + npts_seg-1);
                if strcmp(g.method,'tort')
                    % Tort
                    amplitude_mod(count:count + npts_seg-1) = (g.Ac*((1-mmod(iseg))*sin(2*pi*fm*ti)+1+mmod(iseg))/2).*sin(2*pi*fc*ti) + g.Am*sin(2*pi*fm*ti);%+ noise_amp*rand(1,length(ti));
                else strcmp(g.method,'wiki')
                    amplitude_mod(count:count + npts_seg-1)  = (1+ mmod(iseg).*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod(iseg).*cos(2*pi*fm*ti) + noise_amp.*rand(1,length(ti));
                end
                amplitude_mod(count:count + npts_seg-1) = awgn(amplitude_mod(count:count + npts_seg-1),g.snr,'measured');
                % Modulating signal generation
                phase_signal(count:count + npts_seg-1) = cos(2*pi*fm*ti);
%                 amplitude_mod(count:count + npts_seg-1) = g.Ac*(1+mmod(iseg)*cos(2*pi*fm*ti)).*sin(2*pi*fc*ti) + sin(2*pi*fm*ti)  + noise_amp*rand(1,length(ti)); %
                count = count + npts_seg;
            else
                ti = t(count:end);
                phase_signal(count:count+ length(ti)-1) = cos(2*pi*fm*ti);
                if strcmp(g.method,'tort')
                % Tort
                amplitude_mod(count:count+length(ti)-1) = (g.Ac*((1-mmod(iseg))*sin(2*pi*fm*ti)+1+mmod(iseg))/2).*sin(2*pi*fc*ti) + g.Am*sin(2*pi*fm*ti) + noise_amp*rand(1,length(ti));
                else strcmp(g.method,'wiki')
                    amplitude_mod(count:count+length(ti)-1) = (1+ mmod(iseg).*cos(2*pi*fm*ti)).*g.Ac.*sin(2*pi*fc*ti) + mmod(iseg).*cos(2*pi*fm*ti) + noise_amp.*rand(1,length(ti));
                    
                end
                amplitude_mod(count:count+length(ti)-1) = awgn(amplitude_mod(count:count+length(ti)-1), g.snr,'measured');
%                 amplitude_mod(count:count+length(ti)-1) = g.Ac*(1+mmod(iseg)*cos(2*pi*fm*ti)).*sin(2*pi*fc*ti)+ sin(2*pi*fm*ti) + noise_amp*rand(1,length(ti));
            end
        end
    end
else
    % Modulating signal generation
    phase_signal = cos(2*pi*fm*t);
    % Tort
    m = -m+1;
    amplitude_mod = (g.Ac*((1-mmod)*sin(2*pi*fm*ti)+1+mmod)/2).*sin(2*pi*fc*ti) + g.Am*sin(2*pi*fm*ti) + 10*rand(1,length(t));
    amplitude_mod = awgn(amplitude_mod, g.snr,'measured');
    
%         amplitude_mod = g.Ac*(1+mmod*cos(2*pi*fm*t)).*sin(2*pi*fc*t)+ sin(2*pi*fm*ti) + 10*rand(1,length(t));
end

if g.plot_flag
    figure('Units','normalized','Position',[0.3049 0.2544 0.5486 0.6311]);
    ax(1) = subplot(3,1,1);
    plot(t,phase_signal);
    title (['Modulator/Phase f = ' num2str(fm) 'Hz']);
    %ylabel ('Amplitude');
    ax(1).XTickLabel = [];
    grid on;
    
    
    % Plot 2. Carrier
    ax(2) = subplot(3,1,2);
    plot(t,carrier_signal);
    title (['Carrier/Amplitude f = ' num2str(fc) 'Hz']);
    %ylabel ('Amplitude (au)');
    ax(2).XTickLabel = [];
    grid on;
    
    % Plot 3. Modulated signal
    ax(3) =subplot(3,1,3);
    plot(t,amplitude_mod); hold on;
    plot(t,m*max(amplitude_mod),'color',[215 25 25]./255,'linewidth',2)
    title ('Amplitude Modulated Signal');
    xlabel ('time(sec)');
%     ylabel ('Amplitude (au)');
    grid on;
    hlegend = legend('Signal', 'Modulation');
    hlegend.Box = 'off';
    hlegend. Position = [0.755 0.2650 0.1032 0.0431];
    
    linkaxes(ax,'x');
    
end
