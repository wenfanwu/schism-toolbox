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
%       bathymetry or boundary lines
% disp_schism_hgrid(Mobj, opt_flags, varargin) specifies the visualization styles.
% 
%% Example
% clc;clearvars
% mesh_file  = 'Exp1_BYES\inputs\BYES_64676.mat';
% Mobj = mesh2schism(meshfile);
% figure; disp_schism_hgrid(Mobj, [1 0])
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
% opt_flags - option flags; numeric
%       a two-element vector to determine whether to show the bathymetry
%       and boundary lines. Default: opt_flags = [1 1]; both elements can
%       be specified as 0/1/2. 
% varargin - name-value inputs;
%       the name-value inputs inherited from the 'patch' function.
%
%% Output Arguments
% None
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 20 Aug 2025
% Email: wwu@vims.edu
% 
% See also: mesh2schism and disp_schism_var

%% Parse inputs
if nargin < 2; opt_flags = [1 0]; end
dep_flag = opt_flags(1); bnd_flag = opt_flags(end);

opts_def = {'EdgeColor', 'k','Linewidth',0.0125,'EdgeAlpha', 0.05};
varargin = [opts_def(:)', varargin(:)'];

tri = Mobj.tri;
if numel(find(Mobj.i34==4)) == 0
    tri = tri(:, 1:3);
else
    tri = fillmissing(tri,"previous",2); % for speeding
end
%% Bathymetry
cmap = jet(25);
switch dep_flag
    case 0 % none color
        patch('Faces', tri,'Vertices', [Mobj.lon(:) Mobj.lat(:)], 'FaceColor', 'none', varargin{:});
        
    case 1 % interpolated
        patch('Faces', tri,'Vertices', [Mobj.lon(:) Mobj.lat(:)], 'FaceVertexCData', -Mobj.depth(:), 'FaceColor','interp', varargin{:});
        cbar = colorbar;
        cbar.Label.String = 'depth (m)';
        colormap(cmap)
        dcm = datacursormode;
        dcm.UpdateFcn = @schism_datatips;
        dcm.Enable = 'off';

    case 2 % color on the mesh
        trimesh(tri, Mobj.lon(:), Mobj.lat(:), -Mobj.depth(:));
        view(0, 90)
        cbar = colorbar;
        cbar.Label.String = 'depth (m)';
        colormap(cmap)
        dcm = datacursormode;
        dcm.UpdateFcn = @schism_datatips;
        dcm.Enable = 'off';

end
box on;
xlabel('Longitude', 'FontWeight','bold')
ylabel('Latitude', 'FontWeight','bold')
hold on

%% Boundary
switch bnd_flag
    case 0

    case {1, 2}
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

        if bnd_flag==2 % add boundary index
            for sn = 1:Mobj.obc_counts
                obc_nodes = Mobj.obc_nodes(:,sn);
                obc_nodes(obc_nodes==0) = [];
                sx = Mobj.lon(obc_nodes);
                sy = Mobj.lat(obc_nodes);
                idx = geomin(sx, sy, mean(sx), mean(sy), 1, Mobj.coord);
                text(sx(idx), sy(idx), [' ', num2str(sn), ' '], 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                    'FontWeight', 'bold', 'BackgroundColor', 'w','Margin', 0.001, 'EdgeColor', 'k', 'Color', 'b')
            end
        end
end
axis image
hold off
end