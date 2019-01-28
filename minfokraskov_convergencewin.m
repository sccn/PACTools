function  [I,Ilocal,kconv,difvarvect, totalIlocal] = minfokraskov_convergencewin(Xorig,Yorig,varargin)
% ADD HELP here
% Xorig,Yorig colums vectors for single trials

if nargin < 3
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
try g.k0;                catch, g.k0  = 1;                        end;
try g.k;                 catch, g.k   = [];                       end;
try g.kraskovmethod;     catch, g.kraskovmethod   = 1;            end;
try g.xdistmethod;       catch, g.xdistmethod     = 'seuclidean'; end;
try g.ydistmethod;       catch, g.ydistmethod     = 'seuclidean'; end;
try g.jointdistmethod;   catch, g.jointdistmethod = 'seuclidean'; end;
try g.yvarnorm_circ;     catch, g.yvarnorm_circ   = 0;            end;
try g.xvarnorm_circ;     catch, g.xvarnorm_circ   = 0;            end;
try g.varthresh;         catch, g.varthresh       = 0.05;         end;
try g.kstep;             catch, g.kstep           = 1;            end;
try g.saveItot;          catch, g.saveItot        = 0;            end;
try g.maxk;              catch, g.maxk            = 40;           end;

% Check dimension of signals
X         = Xorig(:);
Y         = Yorig(:);
nlat      = size(X,1);
nlatsorig = size(Xorig,1);

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

dxnorm = dxnorm/max(max(dxnorm));
dynorm = dynorm/max(max(dynorm));
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
totalIlocal = [];

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
    if g.saveItot
         if g.k ==1
             totalIlocal = Ilocal;
         else
             totalIlocal(counter,:)= Ilocal;
         end
             
        %totalIlocal{counter}= Ilocal;
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