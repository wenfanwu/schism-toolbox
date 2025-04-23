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
%       1) bnd_width (km) determines the width of max-nudging zone adjacent to the boundary; 
%       2) cutoff_dist (km) determines the width of total nudging zone at the boundary; 
%       3) nf_max is the maximum nudging factor within the max-nudging zone.
%       bnd_width should be less than cutoff_dist normally. 
%       default: ndg_inputs = [20, 60, 4e-5]; 
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
% To put it simply, the nudging factor will keep nf_max within the
% max-nudging zone, and then decreases sharply away from the open boundary,
% of which the decreasing gradient is controlled by bnd_width/cutoff_dist.,
% with a higher ratio means a stronger decreasing gradient.   
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 29 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: write_schism_nu_nc

%% Parse inputs
if nargin < 2; nudge_inputs = [20, 60, 4e-5]; end
if nargin < 3; obc_bnds = 1:Mobj.obc_counts; end
if nargin < 4; disp_flag = 'on'; end

bnd_width = nudge_inputs(1); 
cutoff_dist = nudge_inputs(2);
nf_max = nudge_inputs(3);
if cutoff_dist<=bnd_width; error('the cutoff_dist must be greater than bnd_width!'); end

%% Calculation
uy = Mobj.lat; ux = Mobj.lon;
obc_nodes = Mobj.obc_nodes(:, sort(obc_bnds));
obc_nodes(obc_nodes==0) = [];
obc_nodes = obc_nodes(:);

if strncmpi(Mobj.coord, 'geographic', 3)
    fcn = @(idx) distance(Mobj.lat(idx), Mobj.lon(idx), Mobj.lat, Mobj.lon, [6378.137 0.0818191910428158])';  % km
    dist_tmp = arrayfun(fcn, obc_nodes, 'UniformOutput', false);
else
    dist_tmp = arrayfun(@(x,y) hypot(ux-x, uy-y)', ux(obc_nodes), uy(obc_nodes), 'UniformOutput', false);
end

dist = cell2mat(dist_tmp); dist = min(dist); 
dist0 = max(dist-bnd_width, 0).^(0.75);  % more smooth 
nudge_factor = (1-tanh(5*dist0./(cutoff_dist-bnd_width))).*nf_max;  % nudge_factor will approach to zero when dist equals to cutoff_dist
nudge_factor(dist>cutoff_dist) = 0;

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