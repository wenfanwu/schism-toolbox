function Mobj = read_schism_vgrid(Mobj, vgrid_file, version_num)
% Read the vertical grids from vgrid.in file.
%
%% Syntax
% Mobj = read_schism_vgrid(Mobj, vgrid_file, version_num)
%
%% Description
% Mobj = read_schism_vgrid(Mobj, vgrid_file, version_num) reads vertical
% grid from an existed vgrid.in file.
%
%% Example
% vgrid_file = 'Exp1_BYS/inputs/vgrid.in';
% Mobj = read_schism_vgrid(Mobj, vgrid_file, 'v5.10')
%
%% Input Arguments
% Mobj --- the mesh object with no vertical layers
% vgrid_file --- the absolute filepath of vgrid.in file.
% version_num --- the version num of SCHISM ('v5.10' or 'v5.9'); This
% option is needed for the LSC2 grid, as the format of vgrid.in has changed
% since v5.10 for LSC2 grid. v5.10 is default.
%
%% Output Arguments
% Mobj --- the mesh object with vertical layers added.
%
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Instiute of Marine Science in 2023.
% Last Updated on 24 Apr 2025
% Email: wwu@vims.edu
%
% See also: read_schism_hgrid

%% Parse inputs
if nargin < 3; version_num = 'v5.10'; end

%% Load vgrid.in file
fid = fopen(vgrid_file, 'r'); 
vdata = textscan(fid, '%s', 'Delimiter', '\n'); vdata = vdata{1}; 
fclose(fid);

% determine the vgrid type (LSC2 or SZ)
ind_exc = find(vdata{1} == '!', 1) - 1;
if isempty(ind_exc)
    ivcor = str2num(vdata{1});  %#ok<*ST2NM>
else
    ivcor = str2num(vdata{1}(1:ind_exc));
end

% read vertical grids
switch ivcor
    case 1
        ind_exc = find(vdata{2} == '!', 1) - 1;
        if isempty(ind_exc)
            maxLev = str2num(vdata{2});
        else
            maxLev = str2num(vdata{2}(1:ind_exc));
        end
        vdata(1:2) = [];

        if strcmpi(version_num, 'v5.9') % old format of vgrid.in
            disp('read vgird.in (LSC2) for the version v5.9 and below')
            nNodes = numel(vdata);
            vgrids = nan(maxLev, nNodes);
            for ii = 1:nNodes
                vals = sscanf(vdata{ii}, '%f');
                zvals = sort(vals(3:end), 'descend');
                if zvals(1) ~= 0 || zvals(end) ~= -1
                    error('failed to read the vgrid.in due to split error!')
                end
                vgrids(1:numel(zvals), ii) = zvals;
            end
        end

        if strcmpi(version_num, 'v5.10') % new format of vgrid.in
            disp('read vgird.in (LSC2) for the version v5.10 and above')
            raw_matrix = cell2mat(cellfun(@(x) sscanf(x, '%f')', vdata(2:maxLev+1), 'UniformOutput', false));
            vgrids = flipud(raw_matrix(:, 2:end));
            vgrids(vgrids == -9) = nan;
        end
        
        nLevs = sum(~isnan(vgrids)); 
        maxLev = max(nLevs);
        vgrids = vgrids(1:maxLev, :);    % trim unused rows
        depLayers = Mobj.depth' .* vgrids;  % auto-broadcasting (R2016b+)
        disp(['the mean # of levels is ', num2str(mean(nLevs), '%.2f')])

        Mobj.nLevs = nLevs;
        Mobj.maxLev = maxLev;
        Mobj.vgrids = vgrids;
        Mobj.depLayers = depLayers;
        Mobj.vtype = 'LSC2';
    case 2
        disp('read vgird.in (SZ)')
        warning on; warning('not work yet')
end

end












