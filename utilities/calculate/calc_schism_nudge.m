function [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs, obc_bnds, disp_flag)
% Calculate open boundary nudging factors
%
%% Syntax
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs)
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs, obc_bnds)
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs, obc_bnds, disp_flag)
%
%% Description
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs)
%       specifies the paramters controlling nudging factors
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs, obc_bnds)
%       specifies the open boundaries to be nudging 
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj, nudge_inputs, obc_bnds, disp_flag)
%       determines whether to show the results or not.
%
%% Examples 
% [nudge_factor, nudge_nodes] = calc_schism_nudge(Mobj,  [20, 60, 4e-5], 1:Mobj.obc_counts, 'on');
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing mesh info.
% nudge_inputs - nudge inputs; numeric
%       three parameters controlling the nudge factors at the open
%       boundary. ngd_inputs = [bnd_width, cutoff_dist, nf_max]; 
%       1) edge_wdith (km) determines the width of max-nudging zone adjacent to the boundary; 
%       2) total_width (km) determines the total width of nudging zone at the boundary; 
%       3) nf_max is the maximum nudging factor within the max-nudging zone.
%       edge_wdith should be less than total_width normally. Default: ndg_inputs = [20, 60, 4e-5]; 
% obc_bnds - open boundaries; numeric
%       the index of open boudaries to be nudged, obc_bnds=1 means that
%       only the first open boundary will be used for nudging. Default:
%       obc_bnds = 1:Mobj.obc_counts.
% disp_flag - display flag; char
%       the flag determining whethe to show the results or not (on/off).
%       default: disp_flag = 'on'.
% 
%% Output Arguments
% nudge_factor - nudging values; double
%       the nudging values@nodes
% nudge_nodes - nudging nodes; double
%       the nodes to be nudging, namely "map_to_global_node".
% 
%% Notes
% The nudging factor stays at nf_max within bnd_width from the boundary,
% then decreases linearly. The decay rate is controlled by bnd_width/cutoff_dist:
% a larger ratio results in a steeper drop.
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 4 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_nu_nc

%% Parse inputs
if nargin < 2; nudge_inputs = [20, 60, 4e-5]; end
if nargin < 3; obc_bnds = 1:Mobj.obc_counts; end
if nargin < 4; disp_flag = 'on'; end

edge_width = nudge_inputs(1);  % km 
total_width = nudge_inputs(2);  % km
nf_max = nudge_inputs(3);

if total_width <= edge_width; error('the cutoff_dist must be greater than bnd_width!'); end
%% Calculation
obc_nodes = Mobj.obc_nodes(:, sort(obc_bnds)); 
obc_nodes(obc_nodes==0) = []; obc_nodes = obc_nodes(:);

lon_obc = Mobj.lon(obc_nodes); 
lat_obc = Mobj.lat(obc_nodes);
[lon1, lon2] = meshgrid(lon_obc, Mobj.lon);
[lat1, lat2] = meshgrid(lat_obc, Mobj.lat);

% Distance matrix from all points to all open boundary points.
if strncmpi(Mobj.coord, 'geographic', 3)
    dist_m = haversine_dist(lat1, lon1, lat2, lon2); % Haversine formula (m)
else
    dist_m = hypot(lon1-lon2, lat1-lat2);
end

dist_pts = min(dist_m, [], 2); 
nudge_factor = nf_max*(dist_pts - total_width)/(edge_width - total_width); % linear transition
nudge_factor = min(nf_max, max(nudge_factor, 0)); 
nudge_nodes = find(nudge_factor~=0);

%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, nudge_factor)
    hold on
    plot_schism_bnds(Mobj)
    axis image
    scatter(Mobj.lon(nudge_nodes), Mobj.lat(nudge_nodes), 3, 'o', 'filled', 'k')
    colormap(jet(25))
    auto_center
end

end
function hdist = haversine_dist(lat1, lon1, lat2, lon2)
% Haversine distance with standard sphere radius.

R = 6378.137; % earch radius (km)
dlat = deg2rad(lat2 - lat1); dlon = deg2rad(lon2 - lon1);
a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
hdist = R .* c; % km

end