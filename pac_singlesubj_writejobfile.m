% Author:  Ramon Martinez-Cancino, UCSD, INC , SCCN  2019
%
% This project was in part supported by the European Union's Horizon 2020
% research and innovation programme under Marie Sklodowska-Curie Action
% agreement no. 750947 (BIREHAB)
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

tmpopt = options(5:end);
optname = tmpopt(1:2:end);
optvals =  tmpopt(2:2:end);

opttext = '{';
for i = 1:length(optname) % Firts eight arguments are passed directly to pop_pac
    if isstr(optvals{i})
        opttext = [opttext ' ''''' optname{i} '''''' ' ''''' optvals{i} ''''''];
    elseif isnumeric(optvals{i})
        opttext = [opttext ' ''''' optname{i} '''''' ' [' num2str(optvals{i}) ']'];
    elseif islogical(optvals{i})
        opttext = [opttext ' ''''' optname{i} '''''' ' ' num2str(optvals{i})];
    end
end
opttext = [opttext '}'];

fid = fopen( fullfile(tmpJobPath,'pacssnsg_job.m'), 'w');
fprintf(fid, 'eeglab;\n');
fprintf(fid, 'EEG  = pop_loadset(''%s'');\n', EEG.filename);
% Defininng variables
fprintf(fid, 'freqs1                  = [%u %u];\n', freqs1);
fprintf(fid, 'freqs2                  = [%u %u];\n', freqs2);
fprintf(fid, 'indexfreqs1             = %u;\n', indexfreqs1);
fprintf(fid, 'indexfreqs2             = %u;\n', indexfreqs2);
fprintf(fid, 'pooldata                = ''%s'';\n', pooldata);
fprintf(fid, 'pacopt     = eval(''%s'');\n', opttext);

% Calling pop_pac
fprintf(fid,'[EEG,com] = pop_pac(EEG,pooldata,freqs1,freqs2,indexfreqs1,indexfreqs2,''compflag'', ''local'',pacopt{:});\n');
fprintf(fid,'pop_saveset(EEG,''filename'', EEG.filename, ''filepath'', pwd);');
fclose(fid);