function write_schism_ic(Mobj, prefix_name, varData)
% Write the *ic file for SCHISM (hvar or vvar)
%
%% Syntax
% write_schism_ic(Mobj, prefix_name, varData)
%
%% Description
% write_schism_ic(Mobj, prefix_name, varData) writes the [prefix_name].ic
%       file using 'varData'.
%
%% Example
% ssh = rand(Mobj.nNodes ,1);
% write_schism_ic(Mobj, 'elev', ssh)  % hvar-type
% write_schism_ic(Mobj, 'salt', 30)  
% ts = [-50:2:0; 25*ones(1,26); 2*ones(1,26)]';
% write_schism_ic(Mobj, 'ts', ts)  % vvar-type
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       prefix name of the *.ic file. e.g. prefix_name = 'elev', or 'elev.ic'.
%       the file extension will be omitted automatically. Both hvar- and
%       vvar-type ic files are supported.
% varData - input data; numeric
%       the initial data to be written. Supported formats:
%       1) scalar or vector (nNodes×1): creates an hvar-type ic file 
%          (e.g., temp.ic, salt.ic, or [mod]_hvar_#.ic). A scalar implies
%          spatially-uniform initial data.
%       2) 2-D matrix (M×N, N>1): creates a vvar-type ic file 
%          (e.g., ts.ic or [mod]_vvar_#.ic), where M represents z layers.
%          The first column indicates the vertical layers.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 19 Apr 2025.
% Email: wwu@vims.edu
%
% See also: write_schism_gr3

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.ic$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '.ic']);

%% Check inputs
if isscalar(varData) || (isvector(varData) && length(varData)==Mobj.nNodes)
    ftype = 'hvar';
    varData = varData(:).*ones(Mobj.nNodes, 1);  % nNodes*1
elseif ismatrix(varData) && min(size(varData))>=2
    ftype = 'vvar';
    if size(varData,1)~=length(varData); varData = varData'; end  % nLevs*nVars
else
    error('Not supported data format!')
end
if any(isnan(varData(:))); error('NaNs were found in the input data!'); end

warning on
switch prefix_name
    case 'ts'
        if size(varData,2)~=3; error('No enough data to create ts.ic!'); end
        if min(varData(:,2))<-2 || max(varData(:,2))>40
            warning('Invalid temp. values have been removed!')
            varData(:,2) = min(40, max(-2, varData(:,2)));
        end
        if min(varData(:,3))<-2 || max(varData(:,3))>40
            warning('Invalid temp. values have been removed!')
            varData(:,3) = min(40, max(-2, varData(:,3)));
        end
    case 'temp'
        if min(varData)<-2 || max(varData)>40
            warning('Invalid temp. values have been removed!')
            varData = min(40, max(-2, varData));
        end
    case 'salt'
        if min(varData)<0 || max(varData)>42
            warning('Invalid salt. values have been removed!')
            varData = min(42, max(0, varData));
        end
end
%% Begin to write
switch ftype
    case 'hvar'  % temp.ic, salt.ic, [mod]_hvar_#.ic
        write_hvar_ic(filepath, varData, Mobj)
    case 'vvar'  % ts.ic, [mod]_vvar_#.ic
        write_vvar_ic(filepath, varData)
end

disp([prefix_name, '.ic has been created successfully!'])
end

function write_hvar_ic(filepath, varData, Mobj)
% Write hvar-type ic file (e.g., temp.ic, salt.ic, [mod]_hvar_#.ic)

head_line = datestr(now, 'mmm/dd/yyyy HH:MM:SS'); %#ok<TNOW1,DATST>

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
end

function write_vvar_ic(filepath, varData)
% Write vvar-type ic file (e.g., ts.ic, [mod]_vvar_#.ic)

[nLevs, nVars] = size(varData);
[~, idx] = sort(varData(:,1), 'ascend');  % starting from bottom
varData = varData(idx,:);

fid = fopen(filepath,'w');
fprintf(fid, '%d  !# of vertical levels\n', nLevs);
fprintf(fid, ['%d', repmat(' %.4f', [1, nVars]),' !z-coord (starting from bottom)\n'], [1, varData(1,:)]');
fprintf(fid, ['%d', repmat(' %.4f', [1, nVars]),'\n'], [(2:nLevs)', varData(2:end,:)]');
fprintf(fid, '\n\n!Notes: (1) code will extrapolate below bottom/above surface; \n');
fprintf(fid, '!       (2) since cubic spline is used for intepolation, make sure the \n');
fprintf(fid, '!           levels are resolved, i.e., no large jumps in z-coord, to \n');
fprintf(fid, '!           avoid non-monotone behavior in between levels. \n');
fclose(fid);

end