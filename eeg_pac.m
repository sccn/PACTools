% eeg_pac() - compute phase-amplitude coupling (power of first input
%         correlation with phase of second). There is no graphical output
%         to this function.
%
% Usage:
%   >> eeg_pac(x,y,srate);
%   >> [crossfcoh, timesout1, freqs1, freqs2, crossfcohall, alltfX, alltfY,crossfcoh_pval, pacstruct] ...
%                     = eeg_pac(x,y,srate,'key1', 'val1', 'key2', val2' ...);
% Inputs:
%    x       = [float array] 2-D data array of size (times,trials) or
%              3-D (1,times,trials)
%    y       = [float array] 2-D or 3-d data array
%    srate   = data sampling rate (Hz)
%
%    Most important optional inputs
%       'alpha'     = Significance level of the statistical test. If
%                     empty no statistical test is done. Empty by Default.                     
%       'bonfcorr'  = Logical. Apply Bonferroni correction to the alpha value
%                     if true. Default [false]
%       'freqs'     = [min max] frequency limits. Default [minfreq 50], 
%                     minfreq being determined by the number of data points, 
%                     cycles and sampling frequency. Use 0 for minimum frequency
%                     to compute default minfreq. You may also enter an 
%                     array of frequencies for the spectral decomposition
%                     (for FFT, closest computed frequency will be returned; use
%                     'padratio' to change FFT freq. resolution).
%       'freqs2'    = [float array] array of frequencies for the second
%                     argument. 'freqs' is used for the first argument. 
%                     By default it is the same as 'freqs'.
%       'method'    = ['mod'|'corrsin'|'corrcos'|'latphase'] modulation
%                     method or correlation of amplitude with sine or cosine of 
%                     angle (see ref). 'laphase' compute the phase
%                     histogram at a specific time and requires the
%                     'powerlat' option to be set.
%       'methodpac' = {'mvlmi', 'klmi', 'glm'} Method to be use
%                     to compute the phase amplitude coupling. 
%                     mvlmi : Mean Vector Length Modulation Index (Canolty et al. 2006)
%                     klmi  : Kullback-Leibler Modulation Index (Tort et al. 2010)
%                     glm   : Generalized Linear Model (Penny et al. 2008)
%                     Default {'glm'}
%       'nbinskl'   = Number of bins to use for the Kullback Leibler
%                     Modulation Index. Default [18].
%       'nboot'     = Number of surrogate data to use. Default [200]
%       'ntimesout' = Number of output times (int<frames-winframes). Enter a 
%                     negative value [-S] to subsample original time by S.
%       'timesout'  = Enter an array to obtain spectral decomposition at 
%                     specific time values (note: algorithm find closest time 
%                     point in dafreqs1ta and this might result in an unevenly spaced
%                     time array). Overwrite 'ntimesout'. {def: automatic}
%       'powerlat'  = [float] latency in ms at which to compute phase
%                     histogram
%       'ptspercent'= Size in percentage of the segments to shuffle 
%                     when creating surrogate data. Default [0.05]
%       'tlimits'   = [min max] time limits in ms. Default [0 number of
%                     time points / sampling rate]
%
%    Optional Detrending:
%       'detrend'   = ['on'|'off'], Linearly detrend each data epoch   {'off'}
%       'rmerp'     = ['on'|'off'], Remove epoch mean from data epochs {'off'}
%
%    Optional FFT/DFT Parameters:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                     If cycles >0: *longest* window length to use. This
%                     determines the lowest output frequency. Note that this
%                     parameter is overwritten if the minimum frequency has been set
%                     manually and requires a longer time window {~frames/8}
%       'padratio'  = FFT-length/winframes (2^k)                    {2}
%                     Multiplies the number of output frequencies by dividing
%                     their spacing (standard FFT padding). When cycles~=0, 
%                     frequency spacing is divided by padratio.
%       'nfreqs1'    = number of output frequencies for modulating signal (phase). 
%                     For FFT, closest computed frequency will be returned. 
%                     Overwrite 'padratio' effects for wavelets.
%                     Default: usfreqs1e 'padratio'.
%       'nfreqs2'    = number of output frequencies for modulated signal (amplitude). 
%                     For FFT, closest computed frequency will be returned. 
%                     Overwrite 'padratio' effects for wavelets.
%                     Default: usfreqs1e 'padratio'.
%       'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                     Note that for obtaining 'log' spaced freqs using FFT, 
%                     closest correspondant frequencies in the 'linear' space 
%                     are returned.
%       'subitc'    = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence 
%                    (ITC) from x and y. This computes the  'intrinsic' coherence
%                     x and y not arising from common synchronization to 
%                     experimental events. See notes. {default: 'off'}
%       'itctype'   = ['coher'|'phasecoher'] For use with 'subitc', see timef()
%                     for more details {default: 'phasecoher'}.
%       'subwin'    = [min max] sub time window in ms (this windowing is
%                     performed after the spectral decomposition).
%       'alltfXstr' = Structure with the TF decomposition. The strcucture
%                     have the fields {alltfs, freqs, timesout}. This is intended to be
%                     used from pop_pac when several channels/components are comuted at a
%                     time. IN this way the TF decompositio does not have to be computed every time.
%       'alltfYstr' = same as 'alltfXstr'
%                     
% Outputs: 
%        crossfcoh      = Matrix (nfreqs1,nfreqs2,length(timesout1)) of phase amplitude 
%                       coupling values.
%        timesout1      = Vector of output times (window centers) (ms).
%        freqs1         = Vector of frequency bin centers for first argument (Hz).
%        freqs2         = Vector of frequency bin centers for second argument (Hz).
%        alltfX         = spectral decomposition of X
%        alltfY         = spectral decomposition of Y
%        crossfcoh_pval = Matrix (nfreqs1,nfreqs2,length(timesout1)) of Pvalues for the
%                        phase amplitude coupling values.
%        pacstruct      = structure containing the paramters and results of the computation             
% Author: Arnaud Delorme, SCCN/INC, UCSD 2005-
%
% Ref: Testing for Nested Oscilations (2008) J Neuro Methods 174(1):50-61
%
% See also: timefreq(), crossf()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [crossfcoh, timesout1, freqs1, freqs2, crossfcohall, alltfX, alltfY,crossfcoh_pval, pacstruct] = eeg_pac(X, Y, srate, varargin);
    
if nargin < 1
    help pac; 
    return; 
end;

crossfcoh_pval = []; pacstruct = []; crossfcohall = [];

% deal with 3-D inputs
% --------------------
if ndims(X) == 3, X = reshape(X, size(X,2), size(X,3)); end;
if ndims(Y) == 3, Y = reshape(Y, size(Y,2), size(Y,3)); end;

frame = size(X,1);
pacmethods_list = {'plv','mvlmi','klmi','glm'} ;

g = finputcheck(varargin, ...
                { 'alpha'         'real'     [0 0.2]                   [];
                  'alltfXstr'     'struct'   struct                    struct;
                  'alltfYstr'     'struct'   struct                    struct;
                  'ptspercent'    'float'    [0 1]                     0.05;
                  'nboot'         'real'     [0 Inf]                   200;
                  'baseboot'      'float'    []                        0;                       %
                  'boottype'      'string'   {'times','trials','timestrials'}  'timestrials';   %
                  'bonfcorr'      'string'   [0 1]                     0;
                  'detrend'       'string'   {'on','off'}              'off';
                  'freqs'         'real'     [0 Inf]                   [0 srate/2];
                  'freqs2'        'real'     [0 Inf]                   [];
                  'freqscale'     'string'   { 'linear','log' }        'linear';
                  'itctype'       'string'   {'phasecoher','phasecoher2','coher'}  'phasecoher';
                  'nfreqs1'       'integer'  [0 Inf]                   [];
                  'nfreqs2'       'integer'  [0 Inf]                   [];
                  'lowmem'        'string'   {'on','off'}              'off';                   %
                  'method'        'string'   { 'mod','corrsin','corrcos'}         'mod';
                  'naccu'         'integer'  [1 Inf]                   250;                     %
                  'newfig'        'string'   {'on','off'}              'on';                    %
                  'padratio'      'integer'  [1 Inf]                   2;
                  'rmerp'         'string'   {'on','off'}              'off';
                  'rboot'         'real'     []                        [];                      %
                  'subitc'        'string'   {'on','off'}              'off';
                  'subwin'        'real'     []                        []; ...
                  'gammapowerlim' 'real'     []                        []; ...                  %
                  'powerlim'      'real'     []                        []; ...                  %
                  'powerlat'      'real'     []                        []; ...                  %
                  'gammabase'     'real'     []                        []; ...                  %
                  'timesout'      'real'     []                        []; ...
                  'methodpac'     'string'   pacmethods_list           'glm';
                  'nbinskl'       'integer'  [1 Inf]                   [];
                  'ntimesout'     'integer'  []                        length(X)/2; ...
                  'tlimits'       'real'     []                        [0 frame/srate];
                  'title'         'string'   []                        '';                      %
                  'vert'          {'real','cell'}  []                  [];
                  'cycles'        'real'     [0 Inf]                   [3 0.5];
                  'cycles2'       'real'     [0 Inf]                   [3 0.5];   
                  'verbose'       'string'   {'on','off'}              'off';
                  'winsize'       'integer'  [0 Inf]                   max(pow2(nextpow2(frame)-3),4) }, 'pac');

if isstr(g), error(g); end;

% more defaults
% -------------
if isempty(g.freqs2),   g.freqs2   = g.freqs;   end;

% remove ERP if necessary
% -----------------------
X = squeeze(X);
Y = squeeze(Y);X = squeeze(X);
trials = size(X,2);
if strcmpi(g.rmerp, 'on')
    X = X - repmat(mean(X,2), [1 trials]);
    Y = Y - repmat(mean(Y,2), [1 trials]);
end;

% perform timefreq decomposition
% ------------------------------
if isempty(fieldnames(g.alltfXstr))
     [alltfX, freqs1, timesout1] = timefreq(X, srate, 'ntimesout',  g.ntimesout, 'timesout',  g.timesout,  'winsize',  g.winsize, ...
                                                      'tlimits',    g.tlimits,   'detrend',   g.detrend,   'itctype',  g.itctype, ...
                                                      'subitc',     g.subitc,    'cycles',    g.cycles,    'padratio', g.padratio, ...
                                                      'freqs',      g.freqs,     'freqscale', g.freqscale, 'nfreqs',   g.nfreqs1,...
                                                      'verbose',    g.verbose);
else
    alltfX    = g.alltfXstr.alltf;
    freqs1    = g.alltfXstr.freqs;
    timesout1 = g.alltfXstr.timesout;
end

if isempty(fieldnames(g.alltfYstr))
    [alltfY, freqs2, timesout2] = timefreq(Y, srate, 'ntimesout',  g.ntimesout, 'timesout',  g.timesout,  'winsize',  g.winsize, ...
                                                     'tlimits',    g.tlimits,   'detrend',   g.detrend,   'itctype',  g.itctype, ...
                                                     'subitc',     g.subitc,    'cycles',    g.cycles2,   'padratio', g.padratio, ...
                                                     'freqs',      g.freqs2,    'freqscale', g.freqscale, 'nfreqs',   g.nfreqs2,...
                                                     'verbose',    g.verbose); 
else
    alltfY    = g.alltfYstr.alltf;
    freqs2    = g.alltfYstr.freqs;
    timesout2 = g.alltfYstr.timesout;
end
 
%%
% check time limits
% -----------------
if ~isempty(g.subwin)
    ind1      = find(timesout1 > g.subwin(1) & timesout1 < g.subwin(2));
    ind2      = find(timesout2 > g.subwin(1) & timesout2 < g.subwin(2));
    alltfX    = alltfX(:, ind1, :);
    alltfY    = alltfY(:, ind2, :);
    timesout1 = timesout1(ind1);
    timesout2 = timesout2(ind2);
end;
if length(timesout1) ~= length(timesout2) | any( timesout1 ~= timesout2)
    if strcmpi(g.verbose, 'on')
        disp('Warning: Time points are different for X and Y. Use ''timesout'' to specify common time points');
    end
    [vals ind1 ind2 ] = intersect_bc(timesout1, timesout2);
    if strcmpi(g.verbose, 'on')
        fprintf('Searching for common time points: %d found\n', length(vals));
    end
    if length(vals) < 10, error('Less than 10 common data points'); end;
    timesout1 = vals;
    timesout2 = vals;
    alltfX = alltfX(:, ind1, :);
    alltfY = alltfY(:, ind2, :);
end;
%%
% scan accross frequency and time
% -------------------------------
cohboot =[];
if numel(size(alltfX)) ==2, ti_loopend = 1; strial_flag = 1; else  ti_loopend = length(timesout1); strial_flag = 0; end
betastmp = []; peakangletmp = []; normpactmp = []; signifmasktmp = []; 
bin_averagetmp = []; surrogate_pactmp = []; compositestmp = []; nbinskltmp = [];

% Apply Bonferroni correction
if g.bonfcorr && ~isempty(g.alpha)
    g.alpha = g.alpha / (length(freqs1) * length(freqs2) * ti_loopend); 
end

if ~strcmpi(g.method, 'latphase')
    for find1 = 1:length(freqs1)
        if strcmpi(g.verbose, 'on')
            fprintf('Progress : %.2f %% \n', 100*(find1-1)/length(freqs1));
        end
        for find2 = 1:length(freqs2)  
            % Insert average here
            for ti = 1:ti_loopend
                
                % get data
                % --------
                if strial_flag
                    tmpalltfx = squeeze(alltfX(find1,:))';
                    tmpalltfy = squeeze(alltfY(find2,:))';
                 else
                     tmpalltfx = squeeze(alltfX(find1,ti,:));
                     tmpalltfy = squeeze(alltfY(find2,ti,:));
                 end

                tmpalltfx = angle(tmpalltfx);
                tmpalltfy = abs(  tmpalltfy);

                % Compute pac
                [pactmp,pval,signiftmp,pacstructtmp] = eeg_comppac(tmpalltfx,tmpalltfy,'method', g.methodpac,'alpha', g.alpha,'ptspercent',g.ptspercent,'nboot',g.nboot, 'nbinskl', g.nbinskl) ;
                                    
                 crossfcoh(find1,find2,ti)      = pactmp;
                 if ~isempty(pval),                       crossfcoh_pval(find1,find2,ti)        = pval;                       end;
                 if ~isempty(signiftmp),                  signifmasktmp(find1,find2,ti)         = signiftmp;                  end;
                 if ~isempty(pacstructtmp.peakangle),     peakangletmp(find1,find2,ti)          = pacstructtmp.peakangle;     end;
                 if ~isempty(pacstructtmp.beta),          betastmp(find1,find2,ti,:)            = pacstructtmp.beta;          end;
                 if ~isempty(pacstructtmp.normpac),       normpactmp(find1, find2, ti)          = pacstructtmp.normpac;       end;
                 if ~isempty(pacstructtmp.bin_average),   bin_averagetmp{find1, find2, ti}      = pacstructtmp.bin_average;   end;
                 if ~isempty(pacstructtmp.surrogate_pac), surrogate_pactmp(find1, find2, ti, :) = pacstructtmp.surrogate_pac; end;
                 if ~isempty(pacstructtmp.composites),    compositestmp(find1,find2, ti, :)     = pacstructtmp.composites;    end;
                 if ~isempty(pacstructtmp.nbinskl),       nbinskltmp(find1,find2, ti)           = pacstructtmp.nbinskl;       end;
            end;
        end;
    end;
    
    % Populating pacstruct
    pacstruct.method            = g.methodpac;
    pacstruct.freqs_phase       = freqs1;
    pacstruct.freqs_amp         = freqs2;
    pacstruct.alpha             = g.alpha;
    pacstruct.ptspercent        = g.ptspercent;
    pacstruct.nboots            = g.nboot;
    pacstruct.srate             = srate;
    pacstruct.bonfcorr          = g.bonfcorr;
    
    pacstruct.pacval            = crossfcoh;
    pacstruct.pval              = crossfcoh_pval;
    pacstruct.signifmask        = signifmasktmp;
    pacstruct.surrogate_pac     = surrogate_pactmp;
    
    pacstruct.peakangle         = peakangletmp;
    pacstruct.betas             = betastmp;
    pacstruct.bin_average       = bin_averagetmp;
    pacstruct.composites        = compositestmp;
    
    pacstruct.nbinskl           = nbinskltmp;
    pacstruct.timesout          = timesout1;
    
elseif strcmpi(g.method, 'latphase')
     fprintf(2,'eeg_pac message: ''latphase'' option is under consruction \n');
%     % this option computes power at a given latency
%     % then computes the same as above (vectors)
%     
%     %if isempty(g.powerlat)
%     %    error('You need to specify a latency for the ''powerlat'' option');
%     %end;
%         
%     gammapower = mean(10*log10(alltfX(:,:,:).*conj(alltfX)),1); % average all frequencies for power
%     if isempty(g.gammapowerlim),
%         g.gammapowerlim = [ min(gammapower(:)) max(gammapower(:)) ];
%     end;
%     fprintf('Gamma power limits: %3.2f to %3.2f\n', g.gammapowerlim(1), g.gammapowerlim(2)); 
%     power = 10*log10(alltfY(:,:,:).*conj(alltfY));
%     if isempty(g.powerlim)
%         for freq = 1:size(power,1)
%             g.powerlim(freq,:) = [ min(power(freq,:)) max(power(freq,:)) ];
%         end;
%     end;
%     for freq = 1:size(power,1)
%         fprintf('Freq %d power limits: %3.2f to %3.2f\n', freqs2(freq), g.powerlim(freq,1), g.powerlim(freq,2)); 
%     end;
%             
%     % power plot
%     %figure; plot(timesout2/1000, (mean(power(9,:,:),3)-mean(power(9,:)))/50);
%     %hold on; plot(linspace(0, length(Y)/srate, length(Y)), mean(Y'), 'g');
% 
%     % phase with power
%     % figure; plot(timesout2/1000, (mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50);
%     % hold on; plot(timesout1/1000, (mean(gammapower,3)-mean(gammapower(:)))/100, 'r');
%     %figure; plot((mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50+j*(mean(gammapower,3)-mean(gammapower(:)))/100, '.');
%     
%     matsize               = 32;
%     matcenter             = (matsize-1)/2+1;
%     matrixfinalgammapower = zeros(size(alltfY,1),size(alltfX,3),matsize,matsize);
%     matrixfinalcount      = zeros(size(alltfY,1),size(alltfX,3),matsize,matsize);    
%     
%     % get power indices
%     if isempty(g.gammabase)
%         g.gammabase = mean(gammapower(:));
%     end;
%     fprintf('Gamma power average: %3.2f\n', g.gammabase); 
%     gammapoweradd  = gammapower-g.gammabase;
%     gammapower     = floor((gammapower-g.gammapowerlim(1))/(g.gammapowerlim(2)-g.gammapowerlim(1))*(matsize-2))+1;
%     phaseangle     = angle(alltfY);
%     posx           = zeros(size(power));
%     posy           = zeros(size(power));
%     for freq = 1:length(freqs2)
%         fprintf('Processing frequency %3.2f\n', freqs2(freq));
%         power(freq,:,:) = (power(freq,:,:)-g.powerlim(freq,1))/(g.powerlim(freq,2)-g.powerlim(freq,1))*(matsize-3)/2+1;
%         complexval      = power(freq,:,:).*exp(j*phaseangle(freq,:,:));
%         posx(freq,:,:)  = round(real(complexval)+matcenter);
%         posy(freq,:,:)  = round(imag(complexval)+matcenter);
%         for trial = 1:size(alltfX,3) % scan trials
%             for time = 1:size(alltfX,2)
%                 %matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial)) = ...
%                 %    matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial))+1;
%                 matrixfinalgammapower(freq,trial,posx(freq,time,trial),posy(freq,time,trial)) = ...
%                     matrixfinalgammapower(freq,trial,posx(freq,time,trial),posy(freq,time,trial))+gammapoweradd(1,time,trial);
%                 matrixfinalcount(freq,trial,posx(freq,time,trial),posy(freq,time,trial)) = ...
%                     matrixfinalcount(freq,trial,posx(freq,time,trial),posy(freq,time,trial))+1;
%             end;
%         end;
%         %matrixfinal(freq,:,:,:) = convn(squeeze(matrixfinal(freq,:,:,:)), gs, 'same');
%         %tmpmat = posx(index,:)+(posy(index,:)-1)*64+(gammapower(:)-1)*64*64;
%         matrixfinalcount(freq, find(matrixfinalcount(freq,:) == 0)) = 1;
%         matrixfinalgammapower(freq,:,:,:) = matrixfinalgammapower(freq,:,:,:)./matrixfinalcount(freq,:,:,:);
%     end;
%     
%     % average and smooth
%     matrixfinalgammapowermean = squeeze(mean(matrixfinalgammapower,2));
%     for freq = 1:length(freqs2)
%         matrixfinalgammapowermean(freq,:,:) = conv2(squeeze(matrixfinalgammapowermean(freq,:,:)), gauss2d(5,5), 'same');
%     end;
%     %matrixfinalgammapower = matrixfinalgammapower/size(alltfX,3)/size(alltfX,2);
%     
%     %vect = linspace(-pi,pi,50);    
%     %for f = 1:length(freqs2)
%     %    crossfcoh(f,:) = hist(tmpalltfy(f,:), vect);
%     %end;
%     
%     % smoothing of output image
%     % -------------------------
%     %gs = gauss2d(6, 6, 6);
%     %crossfcoh = convn(crossfcoh, gs, 'same');
%     %freqs1    = freqs2;
%     %timesout1 = linspace(-180, 180, size(crossfcoh,2));
% 
%     crossfcoh    = matrixfinalgammapowermean;
%     crossfcohall = matrixfinalgammapower; 
end;
if ~exist('crossfcohall', 'var')
    crossfcohall = [];
end    