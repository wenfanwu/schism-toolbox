function write_schism_ic(Mobj, prefix_name, varData)
% Write the *ic files for SCHISM
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
%
% write_schism_ic(Mobj, 'salt', 30)
%
%% Input Arguments
% Mobj --- the mesh object
% prefix_name --- prefix name for the ic file.
% varData --- variable vector or scalar. If the input is a scalar, it will
% be expanded to a Mobj.nNodes*1 vector automatically.
%
%% Output Arguments
% None
%
%% Notes
% the *ic files have the same format with *gr3 files.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 6 Nov. 2024. 
% Email: wwu@vims.edu
% 
% See also:write_schism_gr3 and write_schism_hotstart

%% Parse inputs
if isscalar(varData)
    varData = varData*ones(1, Mobj.nNodes);
end
head_line = datestr(Mobj.time(1), 'mmm/dd/yyyy'); %#ok<*DATST>
prefix_name = regexprep(prefix_name, '\.ic$', '');

%% Check
switch prefix_name
    case 'temp'
        if min(varData)<-2 || max(varData)>40
            warning on
            warning('Invalid temp. values have been removed!')
            varData = max(-2, varData);
            varData = min(40, varData);
        end
    case 'salt'
        if min(varData)<-2 || max(varData)>40
            warning on
            warning('Invalid salt. values have been removed!')
            varData = max(9, varData);
            varData = min(42, varData);
        end
end
%% Begin to write
fileName = fullfile(Mobj.aimpath, [prefix_name, '.ic']);
fid = fopen(fileName,'wt');
fprintf(fid, [head_line, '\n']);
fprintf(fid, [num2str(Mobj.nElems),' ',num2str(Mobj.nNodes), '\n']);

% Node Info
node_part = [(1:Mobj.nNodes)', Mobj.lon(:), Mobj.lat(:), varData(:)]';
node_fmt = repmat('%d   %14.6f   %14.6f   %13.7e\n', 1, size(node_part,2));
fprintf(fid, node_fmt, node_part(:));

% Elem Info
elem_part = [(1:Mobj.nElems)', Mobj.i34(:).*ones(Mobj.nElems,1), Mobj.tri]';
elem_part_1d = elem_part(:);
elem_part_1d(isnan(elem_part_1d)) = [];

elem_fmt = repmat('%d %d %d %d %d %d\n', Mobj.nElems, 1);
elem_fmt(Mobj.i34==3, :) = repmat('%d %d %d %d %d   \n', sum(Mobj.i34==3), 1);
elem_fmt = elem_fmt'; elem_fmt_1d = elem_fmt(:);
fprintf(fid, elem_fmt_1d, elem_part_1d);

fclose(fid);

disp([prefix_name, '.ic has been created successfully!'])
end