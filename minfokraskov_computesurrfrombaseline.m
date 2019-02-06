function [surrdata] = minfokraskov_computesurrfrombaseline(Xdim,Xbaseline,Ybaseline,n_surrogate,varargin)
% This function compute mipac for ER epoched data. Here mipac is computed
% using a window of the size of the phase and the significance is computed
% by creating surrogates from the baseline.

if nargin < 5
    help minfokraskov;
    return;
end

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('minfokraskov() error: calling convention {''key'', value, ... } error'); return;
end;

try g.k0;                catch, g.k0              = 1;            end
try g.k;                 catch, g.k               = [];           end
try g.karskovmethod;     catch, g.karskovmethod   = 1;            end
try g.xdistmethod;       catch, g.xdistmethod     = 'seuclidean'; end
try g.ydistmethod;       catch, g.ydistmethod     = 'seuclidean'; end
try g.jointdistmethod;   catch, g.jointdistmethod = 'seuclidean'; end
try g.yvarnorm_circ;     catch, g.yvarnorm_circ   = 0;            end
try g.xvarnorm_circ;     catch, g.xvarnorm_circ   = 0;            end
try g.varthresh;         catch, g.varthresh       = 0.05;         end
try g.kstep;             catch, g.kstep           = 1;            end
try g.maxk;              catch, g.maxk            = 40;           end

try g.filterfreq;        catch, g.filterfreq      = [];           end
try g.srate;             catch, g.srate           = [];           end
try g.butterorder;       catch, g.butterorder     = 6;            end

OptionNames = fieldnames(g);
arg ={};
for i =1:length(OptionNames)
    if ~ismember(OptionNames{i},{'srate','butterorder','filterfreq'})
        arg{end+1} = OptionNames{i};
        arg{end+1} = g.(OptionNames{i});
    end
end

% Generating surrogates
surrdata = zeros([n_surrogate Xdim(1)]); % initialize
parfor si = 1:n_surrogate
    display(['Surrogate... ' num2str(si)]);
    idx1 = 2;
    idx2 = min(Xdim(1)*Xdim(2), length(Xbaseline(:)));
    
    % X (Phase)
    phaserandindx = floor(idx1+(idx2-idx1).*rand(Xdim(1)*Xdim(2),1));
    X_surrogate = reshape(Xbaseline(phaserandindx),Xdim(1),Xdim(2));
    
    % Y (Amp)
    amprandindx = floor(idx1+(idx2-idx1).*rand(Xdim(1)*Xdim(2),1));
    Y_surrogate = reshape(Ybaseline(amprandindx),Xdim(1),Xdim(2));
    
    [Isurrtmp,Ilocal] = minfokraskov_convergencewin(X_surrogate,Y_surrogate,arg{:});
    
    if ~isempty(g.filterfreq)
        [b,a] = butter(g.butterorder,g.filterfreq/(g.srate/2),'low');
        Ilocal = filtfilt(b,a,Ilocal);
    end
    
    surrdata(si,:) = Ilocal';
end
end