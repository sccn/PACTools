% minfokraskov_convergencewin_jidt() - Compute pairwise local mutual information using Kraskov method
%                                 with modifications to achieve convergence-like while
%                                 iterating through values of 'k'. The function also allow 
%                                 the inputt of extended data inorder to
%                                 boost the neighbors count. This function
%                                 use and thus, required the toobox JIDT, developed by J. Lizier.
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
%  'varthresh'     -  [0.01:0.5] Threshold of the percent decrease of variance if
%                     iterative method to compute MI local is used. Default: 0.5               
%  'kstep'         -  [integer] Step to increase the number of neighbors 'k' in case 
%                     iterative method to compute MI local is used. Default: 1                  
%  'maxk'          -  Maximun value of 'k' in case iterative method to compute 
%                     MI local is used. Default 100.
%  'saveallmi'     - [0,1] Flag to output(1) or not(0) all the vectors of
%                    local MI correspoding to each value of 'k'. Default: [0]
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

function  [I,Ilocal,kconv,difvarvect, AllILocal] = minfokraskov_convergencewin_jidt(Xorig,Yorig,varargin)
% ADD HELP here
% Xorig,Yorig colums vectors for single trials

kconv = []; difvarvect= []; AllILocal = [];

if nargin < 3
    help minfokraskov_convergencewin_jidt;
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
    disp('minfokraskov_convergencewin_jidt() error: calling convention {''key'', value, ... } error'); return;
end
try g.k0;                catch, g.k0              = 1;            end
try g.k;                 catch, g.k               = [];           end
try g.kraskovmethod;     catch, g.kraskovmethod   = 1;            end
try g.varthresh;         catch, g.varthresh       = 0.05;         end
try g.kstep;             catch, g.kstep           = 1;            end
try g.saveallmi;         catch, g.saveallmi       = 0;            end
try g.maxk;              catch, g.maxk            = 40;           end


% Initializations
flagk = 1;
if isempty(g.k)
    g.k = g.k0;
    flagk = 0;
end

counter     = 1;
difvar      = Inf;
difvarvect  = [];
AllILocal   = [];

% Check dimension of signals
source         = Xorig(:);
destination    = Yorig(:);

% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1');

while difvar >= g.varthresh && g.k < g.maxk
    if flagk, difvar = 0; end
    % 2. Set any properties to non-default values:
    calc.setProperty('k', num2str(g.k));
    
    % 3. Initialise the calculator for (re-)use:
    calc.initialise();
    
    % 4. Supply the sample data:
    calc.setObservations(source, destination);
    
    % 5. Compute the estimate:
    result = calc.computeLocalOfPreviousObservations();
    
    milocal = reshape(result(:),[size(Xorig,1), size(Xorig,2)]);
    Ilocal = milocal(:,1)';
    I      = mean(Ilocal);
    
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