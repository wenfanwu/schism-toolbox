function Mobj = read_schism_vgrid(Mobj, vgrid_file, version_num)
% Read the vertical grids
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
% Last Updated on 24 Oct 2024
% Email: wwu@vims.edu
%
% See also: read_schism_hgrid

%% Parse inputs
if nargin < 3
    version_num = 'v5.10';
end
%% Extract
vdata = importdata(vgrid_file, '%/s',inf);
ind_exc = min(strfind(vdata{1}, '!'))-1;
if isempty(ind_exc)
    ivcor = str2num(vdata{1});  %#ok<*ST2NM>
else
    ivcor = str2num(vdata{1}(1:ind_exc));
end

switch ivcor
    case 1
        ind_exc = min(strfind(vdata{2}, '!'))-1;
        if isempty(ind_exc)
            maxLev = str2num(vdata{2});
        else
            maxLev = str2num(vdata{2}(1:ind_exc));
        end
        vdata(1:2) = [];

        if strcmpi(version_num, 'v5.9') % old format of vgrid.in
            disp('read vgird.in (LSC2) for the version v5.9 and below')
            nNodes = length(vdata);
            vgrids = nan(maxLev, nNodes);
            for ii = 1:nNodes
                tline = double(split(string(strtrim(vdata{ii}))));
                tline = sort(tline(3:end), 'descend');
                if tline(1)~=0 && tline(end)~=-1
                    error('failed to read the vgrid.in due to split error!')
                end
                vgrids(1:length(tline), ii) = tline;
            end
        end

        if strcmpi(version_num, 'v5.10') % new format of vgrid.in
            disp('read vgird.in (LSC2) for the version v5.10 and above')
            ind_btm = double(split(string(strtrim(vdata{1}))));
            nNodes = numel(ind_btm);
            vgrids = nan(maxLev, nNodes);
            for ii = 2:maxLev+1
                tline = double(split(string(strtrim(vdata{ii}))));
                vgrids(maxLev+2-ii,:) = tline(2:end);
            end
            vgrids(vgrids==-9) = nan;
        end

        Mobj.nLevs = sum(~isnan(vgrids));
        Mobj.maxLev = max(Mobj.nLevs);
        disp(['the mean # of levels is ', num2str(mean(Mobj.nLevs), '%.2f')])

        Mobj.master_grid = 'unknown';
        vgrids = vgrids(1:Mobj.maxLev,:);

        Mobj.depLayers = bsxfun(@times, Mobj.depth', vgrids);
        Mobj.vgrids = vgrids;
    case 2
        disp('read vgird.in (SZ)')
        warning('Not work now')
end

end












