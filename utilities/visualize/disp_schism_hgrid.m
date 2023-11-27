function disp_schism_hgrid(Mobj, opt_flags, varargin)
% Visualize horizontal grids.
% 
%% Syntax
% disp_schism_hgrid(Mobj)
% disp_schism_hgrid(Mobj, opt_flags)
% disp_schism_hgrid(Mobj, opt_flags, varargin)
%
%% Description 
% disp_schism_hgrid(Mobj) plots the horizontal grids
% disp_schism_hgrid(Mobj, opt_flags) determines whether to display the
% bathymetry or boundary lines
% disp_schism_hgrid(Mobj, opt_flags, varargin) specifies the types
% 
%% Example
% clc;clearvars
% mesh_file  = 'Exp1_BYES\inputs\BYES_64676.mat';
% Mobj = mesh2schism(meshfile);
% figure
% disp_schism_hgrid(Mobj, [1 0])
%
%% Input Arguments
% Mobj --- the mesh object
% opt_flags --- a two-element vector to determine whether to show
% the bathymetry and boundary lines. Default: opt_flags = [1 1]; the first
% element can be 0/1/2; while the second input can be 0/1; 
% varargin --- Name-value inputs, inherited from the 'patch' function; this
% option becomes invalid when opt_flags(1) equals to 2;
%
%% Output Arguments
% None
% 
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: mesh2schism and disp_schism_var

%% Parse inputs
if nargin < 2
    opt_flags = [1 0];
end
dep_flag = opt_flags(1);
bnd_flag = opt_flags(end);

opts_def = {'EdgeColor', 'k','Linewidth',0.0125,'EdgeAlpha', 0.05};
varargin = [opts_def(:)', varargin(:)'];

%% Bathymetry
cmap = jet(25);

switch dep_flag
    case 0 % none color
        patch('Faces',Mobj.tri,'Vertices', [Mobj.lon(:) Mobj.lat(:)], 'FaceColor', 'none', varargin{:});
    case 1 % interpolated
        patch('Faces', Mobj.tri,'Vertices', [Mobj.lon(:) Mobj.lat(:)], 'FaceVertexCData', -Mobj.depth(:), 'FaceColor','interp', varargin{:});
        cbar = colorbar;
        cbar.Label.String = 'depth (m)';
        colormap(cmap)
    case 2 % color on the mesh
        trimesh(Mobj.tri, Mobj.lon(:), Mobj.lat(:), -Mobj.depth(:));
        view(0, 90)
        cbar = colorbar;
        cbar.Label.String = 'depth (m)';
        colormap(cmap)
end
box on;
xlabel('Longitude (°E)', 'FontWeight','bold')
ylabel('Latitude (°N)', 'FontWeight','bold')
hold on
%% Boundary
switch bnd_flag
    case 0

    case 1
        ind_tmp = Mobj.obc_nodes_tot;
        lon_obc = Mobj.lon(ind_tmp);
        lat_obc = Mobj.lat(ind_tmp);

        ind_tmp = Mobj.land_nodes_tot;
        lon_land = Mobj.lon(ind_tmp);
        lat_land = Mobj.lat(ind_tmp);

        ind_tmp = Mobj.island_nodes_tot;
        lon_island = Mobj.lon(ind_tmp);
        lat_island = Mobj.lat(ind_tmp);

        scatter(lon_obc, lat_obc, 30,'s','filled','b')
        scatter(lon_land, lat_land, 3,'filled','r')
        scatter(lon_island, lat_island, 3,'filled','g')
end
%% Enable the datatips
% Tap the node points to display the water depth
% if dep_flag ~= 0
%     s =  scatter(Mobj.lon, Mobj.lat, 1, 'filled','k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
%     row = dataTipTextRow('depth', Mobj.depth);
%     s.DataTipTemplate.DataTipRows(end+1) = row;
% end

axis image
hold off
end