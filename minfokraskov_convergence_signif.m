% minfokraskov_convergence_signif() - Compute pairwise local mutual information using Kraskov method
%                                     with modifications to achieve convergence  while
%                                     iterating through values of 'k'. The function also allow 
%                                     the inputt of extended data in order to boost the neighbors count 
% Usage:
%   >>  pac = minfokraskov_convergence_signif(X,Y,500); 
%
% Inputs:
%  Xorig        - Vector of signal X  (Latency x 1)
%  Yorig        - Vector of signal Y  (Latency x 1)
%  srate        - Sampling rate of X and Y
%
% Optional inputs:
%  'k0'            -  [integer] Number of neighbors for in the Kraskov algotithm.
%                     If iterative method to compute local MI local is used the first 
%                     iteration for all values of 'k' will start in 'k0'. 
%  'k'             -  [integer] Number of neighbors for in the Kraskov algotithm.
%                     If no iterative method to compute local MI is used the magnitude 
%                     will be computed using 'k'number of neighbors
%  'kraskovmethod' -  [1,2] Kraskov method to use in the computtaion on local MI. 
%                     Default: [1]           
%  'xdistmethod'   -  {'circ', 'myeucl','seuclid'} Method to compute distances among samples X.
%                     Default: 'seuclid' (standardized euclidean)
%  'ydistmethod'   -  ('circ', 'myeucl','seuclid') Method to compute distances among samples Y. 
%                     Default: 'seuclid' (standardized euclidean)
%  'xvarnorm_circ' -  [0,1] Flag to activate circular normalization of the
%                     norm for X. May be nececary if circular magnitude like
%                     phase is used
%  'yvarnorm_circ' -  [0,1] Flag to activate circular normalization of the
%                     norm for Y. May be nececary if circular magnitude like
%                     phase is used
%  'varthresh'     -  [0.01:0.5] Threshold of the percent decrease of variance if
%                     iterative method to compute MI local is used. Default: 0.5               
%  'kstep'         -  [integer] Step to increase the number of neighbors 'k' in case 
%                     iterative method to compute MI local is used. Default: 1   
%  'saveallmi'     -  [0,1] Flag to output(1) or not(0) all the vectors of
%                     local MI correspoding to each value of 'k'. Default: [0]
%  'maxk'          -  Maximun value of 'k' in case iterative method to compute 
%                     MI local is used. Default 100.
%  'normmethod'    -  ('norm' 'zscore', 'none ') Normalization method to
%                     use for X and Y. Default: 'norm'
%  'scaledistmat'  -  [0,1] Flag to perform (1) or not (0) scaling of
%                     distance matrix. Default: 0
%  'ptspercent'    -  [0.01:0.5] Size in percentage of data of the segments to shuffle 
%                     when creating surrogate data. Default: [0.05] 
%  'nboot'         -  [Integer] Number of surrogates generated for statistical significance analysis
%  'butterorder'   -  [Integer] Order of Butterworth filter. Default [4]
%  'alpha'         -  [Real] Significance threshold. Default: [0.05]
%  'filterfreq'    -  Lowpass filter cutoff frequency Default: [] 
%            
% Outputs:
%   Iloc_origsurr  -  Local Mutual Information vector [1 x length(X)]
%   kconv          -  k value for convergence (iterative case only)
%   Iloc_sigval    -  Binary vector indicating the significant values of
%                     Local MI [1xlength(X)]
%   Iloc_pval      -  Vector of p values.[ 1xlength(X)]
%   difvarvect     -  Decrease of variance for each value of k used in the
%                     iterative proccess.
%   AllILocal      -  All vectors of Local MI for each 'k' value used in
%                     the interations. Empty if 'saveallmi' is 0
%   surrdata       -  Surrogate data. NOT intended for regular use. 
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2019
%
% Copyright (C) 2019  Ramon Martinez-Cancino,INC, SCCN
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

function [Iloc_origsurr, kconv0, Iloc_sigval, Iloc_pval, difvarvect, AllILocal,surrdata] = minfokraskov_convergence_signif(X,Y,srate,varargin)

if nargin < 4
    help minfokraskov_convergence_signif;
    return;
end

Iloc_sigval = []; surrdata= []; Iloc_pval= [];

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else, g= []; end
catch
    disp('minfokraskov_convergence_signif() error: calling convention {''key'', value, ... } error'); return;
end;

try g.k0;                catch, g.k0              = 1;             end
try g.k;                 catch, g.k               = [];            end
try g.karskovmethod;     catch, g.karskovmethod   = 1;             end
try g.xdistmethod;       catch, g.xdistmethod     = 'seuclidean';  end
try g.ydistmethod;       catch, g.ydistmethod     = 'seuclidean';  end
try g.yvarnorm_circ;     catch, g.yvarnorm_circ   = 0;             end
try g.xvarnorm_circ;     catch, g.xvarnorm_circ   = 0;             end
try g.varthresh;         catch, g.varthresh       = 0.05;          end
try g.kstep;             catch, g.kstep           = 1;             end
try g.saveallmi;         catch, g.saveallmi       = 0;             end
try g.maxk;              catch, g.maxk            = 200;            end
try g.normmethod;        catch, g.normmethod      = 'norm';        end 
try g.scaledistmat;      catch, g.scaledistmat    = 0;             end

try g.ptspercent;        catch, g.ptspercent      = 0.05;          end
try g.nboot;             catch, g.nboot           = 200;           end
try g.butterorder;       catch, g.butterorder     = 6;             end
try g.alpha;             catch, g.alpha           = [];            end
try g.filterfreq;        catch, g.filterfreq      = [];            end
try g.usejidt;           catch, g.usejidt         = 0;             end

% Input stuff
OptionNames = fieldnames(g);
arg ={};
for i =1:length(OptionNames)
    if ~ismember(OptionNames{i},{'ptspercent','butterorder','alpha','filterfreq'})
        arg{end+1} = OptionNames{i};
        arg{end+1} = g.(OptionNames{i});
    end
end

if g.usejidt
    [~,Iloc_orig,kconv0, difvarvect, AllILocal] = minfokraskov_convergencewin_jidt(X,Y,arg{:});
else
    [~,Iloc_orig,kconv0,difvarvect, AllILocal] = minfokraskov_convergencewin(X,Y,arg{:});
end

if ~isempty(g.filterfreq)
     [b,a] = butter(g.butterorder,g.filterfreq/(srate/2),'low');
     Iloc_origsurr = filtfilt(b,a,Iloc_orig');
else
    Iloc_origsurr = Iloc_orig';
end

if ~isempty(g.alpha)
    % surrogate analyses for significance testing
    surrdata = zeros([g.nboot length(X)]); % initialize
    pts_seg = ceil(length(X)*g.ptspercent);
    for si = 1:g.nboot
        
        if mod(length(X),pts_seg) == 0
            nsegm = length(X)/pts_seg;
            
            % X (Phase)
            shuffledata1 = reshape(X,pts_seg,nsegm);
            permarray   = randperm(size(shuffledata1,2));
            X_surrogate = shuffledata1(:,permarray);
            X_surrogate = X_surrogate(:)';
            
            % Y (Amp)
            shuffledata2 = reshape(Y,pts_seg,nsegm);
            permarray   = randperm(size(shuffledata2,2));
            Y_surrogate = shuffledata2(:,permarray);
            Y_surrogate = Y_surrogate(:)';
        else
            nsegm = floor(length(X)/pts_seg);
            r     = rem(length(X),pts_seg);
            
            % X (Phase)
            remainder_surr = X(end-r+1:end);
            rand_idx       = round(rand(1)*length(X));
            shuffledata1   = reshape(X(1:end-r),pts_seg,nsegm);
            permarray      = randperm(size(shuffledata1,2));
            X_surrogate    = shuffledata1(:,permarray);
            if rand_idx > length(X)-r
                X_surrogate    = [remainder_surr(length(X)-rand_idx+1:end)' X_surrogate(1:end) remainder_surr(1:length(X)-rand_idx)'];
            else
                X_surrogate    = [X_surrogate(1:rand_idx) remainder_surr' X_surrogate(rand_idx+1:end)];
            end
            X_surrogate    = X_surrogate(:)';
            
            % Y (Amp)
            remainder_surr = Y(end-r+1:end);
            rand_idx       = round(rand(1)*length(Y));
            shuffledata2   = reshape(Y(1:end-r),pts_seg,nsegm);
            permarray      = randperm(size(shuffledata2,2));
            Y_surrogate    = shuffledata2(:,permarray);
            if rand_idx > length(X)-r
                Y_surrogate    = [remainder_surr(length(X)-rand_idx+1:end) Y_surrogate(1:end) remainder_surr(1:length(X)-rand_idx)];
            else
                Y_surrogate    = [Y_surrogate(1:rand_idx) remainder_surr Y_surrogate(rand_idx+1:end)];
            end
            Y_surrogate    = Y_surrogate(:)';
            
        end 
        
        arg{find(cell2mat(cellfun(@(x) strcmp(x,'k'), arg,'UniformOutput',0)))+1} = kconv0; % Fixing value of 'k' to the one from convergence.
        
%         [Isurrtmp,Ilocal,kconv,difvarvect] = minfokraskov_convergencewin(X_surrogate',Y_surrogate',arg{:});
        
        if g.usejidt
            [~,Ilocal] = minfokraskov_convergencewin_jidt(X_surrogate',Y_surrogate',arg{:});
        else
            [~,Ilocal] = minfokraskov_convergencewin(X_surrogate',Y_surrogate',arg{:});
        end

        if ~isempty(g.filterfreq)            
            [b,a] = butter(g.butterorder,g.filterfreq/(srate/2),'low');
            Ilocal = filtfilt(b,a,Ilocal');
        end
        surrdata(si,:) = Ilocal;
    end
    
    Iloc_zscore = (Iloc_origsurr' - mean(surrdata)) ./ std(surrdata);
    Iloc_pval   = 1-normcdf(abs(Iloc_zscore));
    Iloc_sigval = zeros(size(Iloc_pval));
    Iloc_sigval(abs(Iloc_pval) < g.alpha) = 1;    
end
end