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

try pacstr.method;                catch, pacstr.method         = []; end
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

