% minfokraskov_convergencewin() - Compute pairwise local mutual information using Kraskov method
%                                 with modifications to achieve convergence-like while
%                                 iterating through values of 'k'. The function also allow 
%                                 the inputt of extended data inorder to boost the neighbors count 
% Usage:
%   >>  pac = minfokraskov_convergencewin(X,Y);
%
% Inputs:
%  Xorig        - Vector of signal X
%  Yorig        - Vector of signal Y
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
%  'maxk'          -  Maximun value of 'k' in case iterative method to compute 
%                     MI local is used. Default 100.
%  'saveallmi'     - [0,1] Flag to output(1) or not(0) all the vectors of
%                    local MI correspoding to each value of 'k'. Default: [0]
%  'normmethod'    -  ('norm' 'zscore', 'none ') Normalization method to
%                     use for X and Y. Default: 'norm'
%  'scaledistmat'  - [0,1] Flag to perform (1) or not (0) scaling of
%                    distance matrix. Default: 0
% Outputs:
%   I              - Mutual Information
%   Ilocal         - Local Mutual Information
%   kconv          - k value for convergence(iterative case only)
%   difvarvect     - Decrease of variance for each value of k used in the
%                    iterative proccess.
%   AllILocal      - All vectors of Local MI for each 'k' value used in
%                    the interations. Empty if 'saveallmi'  is 0
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

function  [I,Ilocal,kconv,difvarvect, AllILocal] = minfokraskov_convergencewin(Xorig,Yorig,varargin)
% ADD HELP here
% Xorig,Yorig colums vectors for single trials

if nargin < 3
    help minfokraskov_convergencewin;
    return;
end

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else, g= []; end
catch
    disp('minfokraskov_convergencewin() error: calling convention {''key'', value, ... } error'); return;
end
try g.k0;                catch, g.k0              = 1;            end
try g.k;                 catch, g.k               = [];           end
try g.kraskovmethod;     catch, g.kraskovmethod   = 1;            end
try g.xdistmethod;       catch, g.xdistmethod     = 'seuclidean'; end
try g.ydistmethod;       catch, g.ydistmethod     = 'seuclidean'; end
try g.yvarnorm_circ;     catch, g.yvarnorm_circ   = 0;            end
try g.xvarnorm_circ;     catch, g.xvarnorm_circ   = 0;            end
try g.varthresh;         catch, g.varthresh       = 0.05;         end
try g.kstep;             catch, g.kstep           = 1;            end
try g.saveallmi;         catch, g.saveallmi       = 0;            end
try g.maxk;              catch, g.maxk            = 40;           end
try g.normmethod;        catch, g.normmethod      = 'norm';       end 
try g.scaledistmat;      catch, g.scaledistmat    = 1;            end

% Check dimension of signals
X         = Xorig(:);
Y         = Yorig(:);
nlat      = size(X,1);
nlatsorig = size(Xorig,1);

% Normalizing data
% switch g.normmethod
%     case  'norm'
%         X = X./norm(X);
%         Y = Y./norm(Y); 
%     case 'zscore'
%         X = zscore(X);
%         Y = zscore(Y);
% end

% Compute distances
if strncmpi(g.xdistmethod,'circ',4)
    %     dxnorm = abs(squareform(pdist(X,@distfun)));
    tmparraydist = pdist(X,@distfun);
    dxnorm = abs(squareform(tmparraydist));
    if g.xvarnorm_circ
         [s,s0] = circ_std(dxnorm(:));
         dxnorm = dxnorm/s0;
    end
elseif strncmpi(g.xdistmethod,'myeucl',6)
    dxnorm = sqDistance(X',X');
else
    %     dxnorm = squareform(pdist(X,g.xdistmethod));
    tmparraydist = pdist(X,g.xdistmethod);
    dxnorm = squareform(tmparraydist);
end

if strncmpi(g.ydistmethod,'circ',4)
    %     dynorm = squareform(pdist(Y,@distfun));
    tmparraydist = pdist(Y,@distfun);
    dynorm = squareform(tmparraydist);
    if g.yvarnorm_circ
        [s,s0] = circ_std(dynorm(:));
        dynorm = dynorm/s0;
    end
elseif strncmpi(g.ydistmethod,'myeucl',6)
    dynorm = sqDistance(Y',Y');
else
    %     dynorm = squareform(pdist(Y,g.ydistmethod));
    tmparraydist = pdist(Y,g.ydistmethod);
    dynorm = squareform(tmparraydist);
end
if g.scaledistmat
    dxnorm = dxnorm/max(max(dxnorm));
    dynorm = dynorm/max(max(dynorm));
end
dz = max(dxnorm,dynorm);

dxi = dxnorm;
dxi(logical(find(eye(size(dxi))))) = []; 
dxi = reshape(dxi,length(dxnorm),length(dxnorm)-1);
dyi = dynorm;
dyi(logical(find(eye(size(dyi))))) = []; 
dyi = reshape(dyi,length(dynorm),length(dynorm)-1);
dzi = dz;
dzi(logical(find(eye(size(dzi))))) = []; 
dzi = reshape(dzi,length(dz),length(dynorm)-1);
[~, knntmp] = sort(dzi,2);

% Initializations
flagk = 1;
if isempty(g.k)
    g.k = g.k0;
    flagk = 0;
end

counter     = 1;
difvar      = Inf;
difvarvect  = [];
AllILocal = [];

% Reducing the dim to the one in the first trial
knntmp = knntmp(1:nlatsorig,:);
dxi    = dxi(1:nlatsorig,:);
dyi    = dyi(1:nlatsorig,:);

while difvar >= g.varthresh && g.k < g.maxk
    if flagk, difvar = 0; end;
    
    knn = knntmp(:,g.k);
    if g.kraskovmethod == 1
        maxeps = max(diag(dxi(:,knn)),diag(dyi(:,knn)));
        nx  = sum(dxi < repmat(maxeps,1,size(dxi,2)),2);
        ny  = sum(dyi < repmat(maxeps,1,size(dyi,2)),2);
        
        I      = psi(g.k) - mean(mean(psi(nx + 1)) + mean(psi(ny + 1))) + psi(nlat);
        Ilocal = psi(g.k) - (psi(nx + 1) + psi(ny + 1)) + psi(nlat);
        
    elseif g.kraskovmethod == 2
        nx = sum(dxi < repmat(diag(dxi(:,knn)),1,size(dxi,2)), 2);
        ny = sum(dyi < repmat(diag(dyi(:,knn)),1,size(dyi,2)), 2);
        
        I      = psi(g.k) - 1/g.k - mean(mean(psi(nx)) + mean(psi(ny))) + psi(nlat);
        Ilocal = psi(g.k) - 1/g.k - (psi(nx) + psi(ny)) + psi(nlat);
    end
    
    if counter ~=1
        difvar = 100*abs(Ilocalvark - var(Ilocal))/abs(Ilocalvark);
    end
    difvarvect(counter) = difvar;
    Ilocalvark          = var(Ilocal);
    g.k                 = g.k + g.kstep;
    if g.saveallmi
         if g.k ==1
             AllILocal = Ilocal;
         else
             AllILocal(counter,:)= Ilocal;
         end
    end
    counter = counter + 1;
    
end
kconv = g.k-g.kstep;
end
%--- AUX Functions (end of main function) ---
function D = distfun(XI,XJ)
% From Circular Statistics Toolbox for Matlab
% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


% if size(XI,1)~=size(XJ,1) && size(XI,2)~=size(XJ,2) && length(XJ)~=1
%     error('Input dimensions do not match.')
% end

D = angle(exp(1i*XI)./exp(1i*XJ));
end

function D = sqDistance(X, Y)
D = sqrt(abs(bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y)));
end