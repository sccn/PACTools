function [Iloc_origsurr, kconv0, Iloc_sigval, Iloc_pval, dvarvect, totalIlocalif,surrdata] = minfokraskov_convergence_signif(X,Y,n_surrogate,pts_seg,varargin)
% Compute MI based on kraskov method using the variance reduction loop.
% Provide a filtered MI if requested in parameters. Used for InstMIPAC ONLY

if nargin < 5
    help minfokraskov;
    return;
end

Iloc_sigval = []; surrdata= []; Iloc_pval= [];

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

try g.k0;                catch, g.k0  = 1; end
try g.k;                 catch, g.k   = []; end
try g.karskovmethod;     catch, g.karskovmethod   = 1; end;
try g.xdistmethod;       catch, g.xdistmethod     = 'seuclidean';  end;
try g.ydistmethod;       catch, g.ydistmethod     = 'seuclidean';  end;
try g.jointdistmethod;   catch, g.jointdistmethod = 'seuclidean';  end;
try g.yvarnorm_circ;     catch, g.yvarnorm_circ   = 0;             end;
try g.xvarnorm_circ;     catch, g.xvarnorm_circ   = 0;             end;
try g.filterfreq;        catch, g.filterfreq      = [];            end;
try g.srate;             catch, g.srate           = [];            end;
try g.varthresh;         catch, g.varthresh       = 1;             end;
try g.kstep;             catch, g.kstep           = 1;             end;
try g.saveItot;          catch, g.saveItot        = 1;             end;
try g.maxkprop;          catch, g.maxkprop        = 40;            end;
try g.maxkprop;          catch, g.maxkprop        = 40;            end;
try, g.butterorder;      catch, g.butterorder     = 6;             end;
try, g.alpha;            catch, g.alpha           = [];            end;

[I,Iloc_orig,kconv0,dvarvect, totalIlocalif] = minfokraskov_convergencewin(X,Y,'k0',g.k0...
    ,'kraskovmethod',g.karskovmethod...
    ,'xdistmethod',  g.xdistmethod...
    ,'ydistmethod',  g.ydistmethod...
    ,'varthresh',    g.varthresh...
    ,'saveItot',     g.saveItot);


if ~isempty(g.filterfreq)
%     Iloc_origsurr = eegfilt(Iloc_orig', g.srate, [], g.filterfreq);
     g.butterorder = 6;
     [b,a] = butter(g.butterorder,g.filterfreq/(g.srate/2),'low');
     Iloc_origsurr = filtfilt(b,a,Iloc_orig');
else
    Iloc_origsurr = Iloc_orig;
end

if ~isempty(g.alpha)
    % surrogate analyses for significance testing
    surrdata = zeros([n_surrogate length(X)]); % initialize
    for si = 1:n_surrogate
        
        if mod(length(X),pts_seg) == 0
            nsegm = length(X)/pts_seg;
            
            % X (Phase)
            shuffledata1 = reshape(X,pts_seg,nsegm);
            permarray = randperm(size(shuffledata1,2));
            X_surrogate = shuffledata1(:,permarray);
            X_surrogate = X_surrogate(:)';
            
            % Y (Amp)
            shuffledata2 = reshape(Y,pts_seg,nsegm);
            permarray = randperm(size(shuffledata2,2));
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
                Y_surrogate    = [remainder_surr(length(X)-rand_idx+1:end)' Y_surrogate(1:end) remainder_surr(1:length(X)-rand_idx)'];
            else
                Y_surrogate    = [Y_surrogate(1:rand_idx) remainder_surr' Y_surrogate(rand_idx+1:end)];
            end
            Y_surrogate    = Y_surrogate(:)';
            
        end
        
        [Isurrtmp,Ilocal,kconv,difvarvect] = minfokraskov_convergencewin(X_surrogate',Y_surrogate'...
            ,'k',             kconv0...
            ,'kraskovmethod', g.karskovmethod...
            ,'xdistmethod',   g.xdistmethod...
            ,'ydistmethod',   g.ydistmethod...
            ,'xvarnorm_circ', g.xvarnorm_circ...
            ,'yvarnorm_circ', g.yvarnorm_circ...
            );
        if ~isempty(g.filterfreq)
            % Ilocal = eegfilt(Ilocal', g.srate, [], g.filterfreq);
            
            [b,a] = butter(g.butterorder,g.filterfreq/(g.srate/2),'low');
            Ilocal = filtfilt(b,a,Ilocal');
        end
        surrdata(si,:) = Ilocal;
    end
    
    Iloc_zscore         = (Iloc_origsurr' - mean(surrdata)) ./ std(surrdata);
    Iloc_pval           = 1-normcdf(abs(Iloc_zscore));
    
    Iloc_sigval = zeros(size(Iloc_pval));
    Iloc_sigval(abs(Iloc_sigval) < g.alpha) = 1;    
end

end