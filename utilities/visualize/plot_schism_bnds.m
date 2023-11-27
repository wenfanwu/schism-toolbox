function plot_schism_bnds(Mobj, opt_flags, varargin)
% Plot the land/island boundaries for SCHISM
% 
%% Syntax
% plot_schism_bnds(Mobj, opt_flags)
% plot_schism_bnds(Mobj, opt_flags, varargin)
%
%% Description 
% plot_schism_bnds(Mobj, opt_flags) plots the land/island boundaries on the
% map
% plot_schism_bnds(Mobj, opt_flags, varargin) specifies the properties for
% the boundary line
%
%% Examples
% plot_schism_bnds(Mobj, [1 1], 'Color', 'k', 'LineWidth', 1)
%
%% Input Arguments
% Mobj --- the mesh object
% opt_flags --- a two-value vector used to decide whether to plot the land and
% island boundaries or not (1/0). Default: opt_flags = [1 1];
% varargin --- properties of the boundary lines, the same as the function 'plot'
%
%% Output Arguments
% None
% 
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-05-30.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if nargin<2
    opt_flags = [1 1];
end

%% Display
if opt_flags(1) == 1
    ind_nodes = double(Mobj.land_nodes);
    ind_tmp = [ind_nodes; -1*ones(1, size(ind_nodes, 2))];  % support multiple segments of boundaries
    ind_tmp(ind_tmp==0) = [];

    ind_nans = ind_tmp == -1;
    ind_tmp(ind_tmp==-1) = 1;
    lon_tot = Mobj.lon(ind_tmp);
    lat_tot = Mobj.lat(ind_tmp);

    lon_tot(ind_nans) = nan;
    lat_tot(ind_nans) = nan;

    plot(lon_tot, lat_tot, varargin{:})
end

if opt_flags(2) == 1
    ind_nodes = double(Mobj.island_nodes);
    ind_tmp = [ind_nodes; -1*ones(1, size(ind_nodes, 2))];
    ind_tmp(ind_tmp==0) = [];

    ind_nans = ind_tmp == -1;
    ind_tmp(ind_tmp==-1) = 1;
    lon_tot = Mobj.lon(ind_tmp);
    lat_tot = Mobj.lat(ind_tmp);

    lon_tot(ind_nans) = nan;
    lat_tot(ind_nans) = nan;

    plot(lon_tot, lat_tot, varargin{:})
end

end






