function write_schism_ic(Mobj, prefix_name, varData)
% Write the *ic file for SCHISM
% 
%% Syntax
% write_schism_ic(Mobj, prefix_name, varData)
%
%% Description 
% write_schism_ic(Mobj, prefix_name, varData) writes the [prefix_name].ic
% file using 'varData'.
% 
%% Example
% ssh = rand(Mobj.nNodes ,1);
% write_schism_ic(Mobj, 'elev', ssh)
% write_schism_ic(Mobj, 'salt', 30)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       the file prefix name. e.g. prefix_name = 'drag', or 'drag.gr3'. 
%       the file extension will be omitted automatically.
% varData - input data; numeric
%       the input data to be written. If '' is a scalar, it will
%       be expanded to a vector automatically , and the data at all nodes
%       are the same.in this case.
%
%% Output Arguments
% None
%
%% Notes
% the *ic files have the same format with *gr3 files.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 5 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also:write_schism_gr3 and write_schism_hotstart

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.ic$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '.ic']);
head_line = datestr(now, 'mmm/dd/yyyy HH:MM:SS'); %#ok<TNOW1,DATST>

%% Check
if isscalar(varData); varData = varData*ones(Mobj.nNodes, 1); end
if any(isnan(varData(:))); error('NaNs were found in the input data!'); end

switch prefix_name
    case 'temp'
        if min(varData)<-2 || max(varData)>40
            warning on
            warning('Invalid temp. values have been removed!')
            varData = max(-2, varData);
            varData = min(40, varData);
        end
    case 'salt'
        if min(varData)<9 || max(varData)>42
            warning on
            warning('Invalid salt. values have been removed!')
            varData = max(9, varData);
            varData = min(42, varData);
        end
end
%% Begin to write
fid = fopen(filepath,'w');
fprintf(fid, [head_line, '\n']);
fprintf(fid,'%d %d\n', Mobj.nElems, Mobj.nNodes);

% Node and Elem
fprintf(fid, '%d   %14.6f   %14.6f    %13.7e\n', [(1:Mobj.nNodes)', Mobj.lon(:), Mobj.lat(:), varData(:)]');
fprintf(fid, '%d %d %d %d %d %d\n', [(1:Mobj.nElems)', Mobj.i34(:), Mobj.tri]');  % include NaN values
fclose(fid);

% Remove NaNs
fileText = fileread(filepath);
fid = fopen(filepath, 'w');
fwrite(fid, strrep(fileText, 'NaN', '')); % case-sensitive
fclose(fid);

disp([prefix_name, '.ic has been created successfully!'])
end