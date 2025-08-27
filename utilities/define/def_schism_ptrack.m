function [xps, yps] = def_schism_ptrack(Mobj, nRegs, mtype)
% Define the locations for particle tracking.
% 
%% Syntax
% [xps, yps] = def_schism_ptrack(Mobj)
% [xps, yps] = def_schism_ptrack(Mobj, nRegs)
% [xps, yps] = def_schism_ptrack(Mobj, nRegs, mtype)
% 
%% Description 
% [xps, yps] = def_schism_ptrack(Mobj) defines the released
%       region for particle tracking.
% [xps, yps] = def_schism_ptrack(Mobj, nRegs) specifies the # of
%       released regions. 
% [xps, yps] = def_schism_ptrack(Mobj, nRegs, mtype) specifies
%       the defining method.
%
%% Example
% [xps, yps] = def_schism_ptrack(Mobj)
% [xps, yps] = def_schism_ptrack(Mobj, 3, 'polygon)
% 
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store the mesh info.
% nRegs - # of regions; numeric
%       the # of released regions; Default: nRegs = 1;
% mtype - method type; char
%       the defining method (polygon/points) to select released regions.
%       Default: mtype =  'polygon'.
%
%% Output Arguments
% xps - x-axis coordinates; numeric
%       the x-axis coordinates of selected points;
% yps - y-axis coordinates; numeric
%       the y-axis coordinates of selected points;
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 27 Nov 2023.
% Email: wwu@vims.edu
% 
% See also: write_schism_ptrack

%% Parse inputs
if nargin < 2; nRegs = 1; end
if nargin < 3; mtype = 'polygon'; end

%% Define on the map
xps = cell(nRegs,1); yps = cell(nRegs,1);

hold on
for iReg = 1:nRegs
    disp('draw a polygon on the map to select released regions and press ENTER')
    switch lower(mtype)
        case 'polygon'
            geo_handle = drawpolygon;
            x_roi = geo_handle.Position(:,1)';
            y_roi = geo_handle.Position(:,2)';
            is_in = inpolygon(Mobj.lon, Mobj.lat, x_roi, y_roi);
            xtmp = Mobj.lon(is_in);
            ytmp = Mobj.lat(is_in);
        case 'points'
            pause;
            dataTips = findobj(gcf, 'Type', 'datatip');
            if ~isempty(dataTips)
                select_coords = cell2array(arrayfun(@(x) [dataTips(x).X; dataTips(x).Y], 1:length(dataTips), 'UniformOutput', false));
                xtmp = select_coords(:,1);
                ytmp =  select_coords(:,2);
            else
                xtmp = [];
                ytmp = [];
            end
    end
    xps{iReg} = xtmp(:); yps{iReg} = ytmp(:);
    scatter(xtmp, ytmp, 10, 'filled', 'm')
end
xps = cell2mat(xps); yps = cell2mat(yps);

end





















