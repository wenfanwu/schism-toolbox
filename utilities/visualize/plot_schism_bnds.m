function plot_schism_bnds(Mobj, bnd_flags, varargin)
% Plot the land/island/ocean boundaries for SCHISM
% 
%% Syntax
% plot_schism_bnds(Mobj, bnd_flags)
% plot_schism_bnds(Mobj, bnd_flags, varargin)
%
%% Description 
% plot_schism_bnds(Mobj, bnd_flags) plots land/island boundaries on the map
% plot_schism_bnds(Mobj, bnd_flags, varargin) specifies the properties for
%       the boundary line or the projection flags
%
%% Examples
% figure
% disp_schism_var(Mobj, Mobj.depth, 'Projection', 'on', 'EdgeColor', 'k')
% hold on
% plot_schism_bnds(Mobj, [1 1], 'Projection', 'on', 'Color', 'r', 'LineWidth', 2)
% axis image
%
%% Input Arguments
% Mobj - mesh object; datastruct
% bnd_flags - boundary flags; double
%       a three-value vector used to decide whether to plot the land or
%       island or ocean boundaries (1/0); default: bnd_flags = [1 1 0];
% varargin --- attributes of boundary lines, the same as the function
%       'plot', or determine the projection flags ('on'/'off')
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 28 Oct 2024.
% Email: wwu@vims.edu
% 
% See also: disp_schism_var

%% Parse inputs
if nargin<2
    bnd_flags = [1 1 1];
end

opts_def = {'Color', 'k','Linewidth',1};
varargin = [opts_def(:)', varargin(:)'];

ind = find(strncmpi(varargin, 'projection', 4), 1, 'last' ); % only the last one is valid
if isempty(ind)
    proj_flag = 'off';
else
    proj_flag = varargin{ind+1};
    varargin(ind:ind+1) = [];
end

%% Display
hold on
if bnd_flags(1) == 1  % display land nodes or not
    plot_bnd_nodes(Mobj, Mobj.land_nodes, proj_flag, varargin{:})
end

if bnd_flags(2) == 1  % display island nodes or not
    island_nodes = [Mobj.island_nodes; Mobj.island_nodes(1,:)];  % form a loop
    plot_bnd_nodes(Mobj, island_nodes, proj_flag, varargin{:})
end

if bnd_flags(3) == 1  % display island nodes or not
    plot_bnd_nodes(Mobj, Mobj.obc_nodes, proj_flag, varargin{:})
end
end

function plot_bnd_nodes(Mobj, bnd_nodes, proj_flag, varargin)
% plot boundary nodes

ind_nodes = double(bnd_nodes);
ind_tmp = [ind_nodes; -1*ones(1, size(ind_nodes, 2))];  % support multiple segments of boundaries
ind_tmp(ind_tmp==0) = [];

ind_nans = ind_tmp == -1;
ind_tmp(ind_tmp==-1) = 1;

lon_tot = Mobj.lon(ind_tmp);
lat_tot = Mobj.lat(ind_tmp);

lon_tot(ind_nans) = nan;
lat_tot(ind_nans) = nan;

if strcmpi(proj_flag,'on')
    proj = projcrs(3395, 'Authority', 'EPSG');  % EPSG 3395: Mercator projection
    [lon_tot, lat_tot] = projfwd(proj, lat_tot, lon_tot);
end

plot(lon_tot, lat_tot, varargin{:})
end






