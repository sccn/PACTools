% eeg_comppac() - Compute phase amplitude coupling (PAC) with the time 
%                 series of the phase and the amplitude. Generate the 
%                 statistical significance of the PAC with surrogate data.
%
% Usage:
%   >> eeg_comppac(x,y);
%   >> [pacval,pval,significant,pacstr] ...
%                     = eeg_comppac(x,y,'key1', 'val1', 'key2', val2' ...);
% Inputs:
%    x      - Time series of the phase [ntimes] or [ntrials]
%    y      - Time series of the amplitude [ntimes] or [ntrials]
%
%   Optional inputs
%       'alpha'         - Significance level of the statistical test. If
%                         empty no statistical test is done.
%                         Default [0.05]
%       'method'        - {'mvlmi', 'klmi', 'glm', 'plv'} Method to be use
%                         to compute the phase amplitude coupling. 
%                         mvlmi : Mean Vector Length Modulation Index (Canolty et al. 2006)
%                         klmi  : Kullback-Leibler Modulation Index (Tort et al. 2010)
%                         glm   : Generalized Linear Model (Penny et al. 2008)
%                         plv   : Phase-locking value (Lachaux, 1999)
%                         Default {'glm'}
%       'nbinskl'       - Number of bins to use for the Kullback Leibler
%                         Modulation Index. Default [18].
%       'nboot'         - Number of surrogate data to use. Default [200]
%       'normpac'       - Normalize the PAC according to Penny et al.
%                        (2008). Not defined for klmi. [0,1] Default [0]
%       'ptspercent'    - Size in percentage of the segments to shuffle 
%                         when creating surrogate data. Default [0.05]
% Outputs:
%        pacval         - Phase Amplitude Coulping Value
%        pval           - p value of the pacval
%        significant    - Logical. Significance of the data depending on
%                         alpha
%        pacstr         - Structure containing the parameters and
%                         intermediate steps of the method
%     
%
% Author: Ramon Martinez Cancino and Joseph Heng, SCCN/INC, UCSD 2016
%
%
% See also:
% Copyright (C) 2016 Ramon Martinez Cancino and Joseph Heng, UCSD, INC, SCCN
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

function [pacval,pval,significant,pacstr] = eeg_comppac(x,y,method, varargin)
pacstr = []; pacval = []; pval = []; significant = [];

% Checking arguments and assigning default values
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('eeg_comppac() error: calling convention {''key'', value, ... } error'); return;
end;

try g.nbinskl;        catch, g.nbinskl      = 18;         end;
try g.normpac;        catch, g.normpac      = 0;          end;
try g.alpha;          catch, g.alpha        = [];         end;
try g.verbose;        catch, g.verbose      = 0;          end;
try g.ptspercent;     catch, g.ptspercent   = 0.05;       end;
try g.nboot;          catch, g.nboot        = 200;        end;
if ~isequal(size(x),size(y))
    error('eeg_comppac() error: X and Y must have the same dimensions');
end

% Initializing pacstr structure
pacstr = create_pacstr('alpha', g.alpha);

% Compute the pacval depending on the method
switch method
    case 'plv'
        % Phase Locking Value
        pacval = eeg_plv(x,y);
        % Print results
        if g.verbose
          fprintf('Phase Locking Value = %.3f \n', pacval);
        end
        % Put values of interest in pac structure
        % Normalizing pacval following Penny et al. (2008)
        if g.normpac
            pacstr.normpac = asin(2*pacval-1);
        end
        
    case 'mvlmi'
        % Mean Vector Length modulation index
        [pactmp, m_raw, composites] = eeg_mvlmi(x,y);
        
        % Normalize mraw
        [pval, surr_mean, surr_std,surrogate_pac] = eeg_pacstatistics(pactmp,x,y,method,g.nboot, g.nbinskl,g.ptspercent);
        if ~isempty(g.alpha)
           pacstr.significant   = pval<g.alpha;
           pacstr.pval = pval;
        else
            pval = [];
        end
        normlength           = (abs(m_raw)-surr_mean)/surr_std; 
        normphase            = angle(m_raw);        
        pacstr.peakangle     = normphase;
        pacval               = normlength;
        pacstr.composites    = composites;

        if g.normpac
            pacstr.normpac = log(pacval);
        end
        
        % Apply the same transformation to the surrogate_pac data so they
        % are comparable with the pacval
        if ~isempty(g.alpha)
            surrogate_pac = (abs(surrogate_pac)-surr_mean)/surr_std;
        end
        if g.verbose
            fprintf('Mean Vector Length modulation index = %.3f \n', pacval);
            fprintf('Highest modulated amplitude at the modulating phase %.3f \n',pacstr.peakangle);
        end
        
    case 'klmi'
        % Kullback-Leibler modulation index
        [pacval, peakangle, bin_average, nbins] = eeg_klmi(x,y,g.nbinskl,1, g.verbose);
        % Print results
        if g.verbose
            fprintf('Kullback-Leibler modulation index = %.3f \n', pacval);
            fprintf('Highest modulated amplitude at the modulating phase %.3f \n',peakangle);
        end
        % Put values of interest in pac structure
        pacstr.peakangle   = peakangle;
%         pacstr.nbinskl     = g.nbinskl;
        pacstr.bin_average = bin_average;
        pacstr.nbinskl     = nbins;
        
        if nbins < g.nbinskl,
            disp('eeg_comppac() warning: number of bins for KL computation was reduced');
        end
        
    case 'glm'
        % General Linear Model
        [pacval, beta] = eeg_glm(x,y);
        % Print results
        if g.verbose
            fprintf('beta = [ %.3f , %.3f, %.3f] \n', beta(1), beta(2), beta(3));
            fprintf('variance explained = %.3f \n', pacval);
        end
        % Put values of interest in pac structure
        pacstr.beta = beta;
        % Normalizing pacval following Penny et al. (2008)
        if g.normpac
            pacstr.normpac = atanh(sqrt(pacval));
        end
end

pacstr.pacval = pacval; % Put the pac value in the pac structure

% Compute the statistical value
if ~isempty(g.alpha)
    if isempty(pval)
        [pval, ~, ~, surrogate_pac] = eeg_pacstatistics(pacval,x,y,method,g.nboot, pacstr.nbinskl,g.ptspercent, g.verbose);
    end
    if g.verbose
        fprintf('p value = %.3f \n', pval);
    end
    pacstr.pval          = pval;
    pacstr.significant   = pval<g.alpha;
    pacstr.surrogate_pac = surrogate_pac;
end

end

% -------------------------------------------------------------------------
% eeg_plv
% -------------------------------------------------------------------------
function [plv] = eeg_plv(phaseX,amplitudeY)
% Compute the Phase Locking Value (PLV) between two signals

    % Extract the phase from the amplitude
    phaseY = angle(hilbert(amplitudeY));
    
    % Compute the phase locking value
    plv = abs(mean(exp(1j * (phaseX-phaseY))));
end
% eeg_klmi
% -------------------------------------------------------------------------
function [pacval, peakangle, bin_average, nbins] = eeg_klmi(phase, amplitude, nbinskl,modbin_flag, verbose)
% Compute the Kullback-Leibler modulation index
    
    nbins = nbinskl;
    
    flag = true;
    while flag
        bin_size     = 2*pi/nbins;             % Compute bin size
        bin_average  = zeros(1,nbins);         % Initialize the histogramme

        for i=1:nbins                         %Cycle through the bins   
        %   fill the bins with the mean amplitude for the phases that correspond to that bin
            bin_average(i) = mean(amplitude(wrapTo2Pi(phase) > (i-1)*bin_size & wrapTo2Pi(phase) < i*bin_size));
        end

        % Normalize the bins to obtain a probability distribution
        bin_average = bin_average/(sum(bin_average));

        % If some bins are empty, then the KL distance cannot be computed (log of a null value)
         if sum(bin_average>0)<nbins && modbin_flag
            oldnbins = nbins;
            nbins = ceil(nbins/2);
            if verbose,
                fprintf('Too many bins to compute Kullback Leibler MI. Reducing the number of bins from %d to %d \n', oldnbins, nbins);
            end
         else
            flag = false;
         end
    end
    
    % Compute KL distance
    pacval = (sum(bin_average.*log(bin_average/(1/nbins))))/log(nbins);
    
    % Compute the peak phase
    [~,index] = max(bin_average);
    peakangle = index*bin_size - bin_size/2;
end
% eeg_mvlmi
% -------------------------------------------------------------------------
function [pacval, m_raw, z] = eeg_mvlmi(phase, amplitude)
% Compute the Mean Vector Length Modulation Index    
    z      = amplitude.*exp(1j*phase);
    m_raw  = mean(z);
    pacval = abs(m_raw);
end
% eeg_glm
% -------------------------------------------------------------------------
function [pacval, beta] = eeg_glm(phase, amplitude)
% General linear Model PAC method
     X = [cos(phase), sin(phase) ones(size(phase))];                      % Building Matrix of regressors. Note : glmfit adds a column of 1s
    [beta,~, stats] = glmfit(X,amplitude,'normal','constant','off');      % Fit the GLM
    pacval = 1- sum(stats.resid.^2)/sum((amplitude-mean(amplitude)).^2);  % 1-var(stats.resid)/var(amplitude); % Calculate the explained variance
end
 % Surrogate Analysis
% -------------------------------------------------------------------------    
function [p_value, surr_mean, surr_std, surrogate_pac] = eeg_pacstatistics(pacval, tmpalltfx, tmpalltfy, method,nboot, nbinskl,ptspercent, verbose)
% Compute the statistical value of pacval estimating a normal distribution from the surrogate data

    % Building surrogates
    pts_segm = ceil(length(tmpalltfx)*ptspercent);

    % Breaking data in blocks to create surrogates
    countsegm = 1;
    for isegm = 1:ceil(length(tmpalltfx)/pts_segm);
        if isegm ~= ceil(length(tmpalltfx)/pts_segm)
            indx = (isegm-1)*pts_segm+[1:pts_segm];
            shuffledata1{isegm} = tmpalltfx(indx);
            shuffledata2{isegm} = tmpalltfy(indx);
        else
            shuffledata1{isegm} = tmpalltfx(countsegm:end);
            shuffledata2{isegm} = tmpalltfy(countsegm:end);
        end
        countsegm = countsegm + pts_segm;
    end

    for isurr = 1:nboot
        % X (Phase)
        permarray = randperm(size(shuffledata1,2));
        xtmp = shuffledata1(permarray);
        X_surrogate(:,isurr) = cat(1,xtmp{:});

        % Y (Amp)
        permarray = randperm(size(shuffledata2,2));
        ytmp = shuffledata2(:,permarray);
        Y_surrogate(:,isurr) = cat(1,ytmp{:});
    end

    % Initializing surrogate PAC
    nsurrogate    = size(X_surrogate,2);
    surrogate_pac = zeros(1,nsurrogate);

    % Compute the surrogate pac value for each surrogate data depending on the given method
    for i=1:nsurrogate
        switch method
            case 'plv'
              surrogate_pac(i) = eeg_plv(X_surrogate(:,i),Y_surrogate(:,i));
            case 'mvlmi'
              surrogate_pac(i) = eeg_mvlmi(X_surrogate(:,i),Y_surrogate(:,i));
            case 'klmi'
              surrogate_pac(i) = eeg_klmi(X_surrogate(:,i),Y_surrogate(:,i),nbinskl,0, verbose);
            case 'glm'
              surrogate_pac(i) = eeg_glm(X_surrogate(:,i),Y_surrogate(:,i));
        end
    end
                
    % Compute p_value by estimating a normal distribution from the surrogate data
    [surr_mean, surr_std] = normfit(surrogate_pac);          % Estimating a normal distribution
    norm_pacval           = (pacval - surr_mean)/surr_std;   % Normalize the pacval
    p_value               = 1-normcdf(abs(norm_pacval));     % Compute p_value    

end

% Functions
function pacstr = create_pacstr(varargin)
% Create a 'pacstr' structure with the values provided.If no inputs
% provided 'pacstr' is created with empty values ([])

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            pacstr.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('create_pacstr() error: calling convention {''key'', value, ... } error'); return;
end;

try pacstr.pacval;                catch, pacstr.pacval         = []; end
try pacstr.pval;                  catch, pacstr.pval           = []; end
try pacstr.peakangle;             catch, pacstr.peakangle      = []; end
try pacstr.beta;                  catch, pacstr.beta           = []; end
try pacstr.normpac;               catch, pacstr.normpac        = []; end
try pacstr.nbinskl;               catch, pacstr.nbinskl        = []; end
try pacstr.alpha;                 catch, pacstr.alpha          = []; end
try pacstr.significant;           catch, pacstr.significant    = []; end
try pacstr.bin_average;           catch, pacstr.bin_average    = []; end
try pacstr.surrogate_pac;         catch, pacstr.surrogate_pac  = []; end
try pacstr.phaseangles;           catch, pacstr.phaseangles    = []; end
try pacstr.amplitudes;            catch, pacstr.amplitudes     = []; end
try pacstr.composites;            catch, pacstr.composites     = []; end
try pacstr.kconv;                 catch, pacstr.kconv          = []; end
try pacstr.diffvar;               catch, pacstr.diffvar        = []; end

end