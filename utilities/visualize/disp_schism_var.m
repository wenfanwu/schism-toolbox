function disp_schism_var(Mobj, varData, varargin)
% Visualize variables on an unstructured grid.
%
%% Syntax
% disp_schism_var(Mobj, varData)
% disp_schism_var(Mobj, varData, varargin)
%
%% Description
% disp_schism_var(Mobj, varData) visualizes variable on an unstructured grid.
% disp_schism_var(Mobj, varData, varargin) specifies the options or projections.
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
%       a datastruct containing the mesh info.
% varData - variable data; double
%       variable vector defined@nodes or elems or edges.
% varargin - input options
%       this function shares the same options with the built-in fuction
%       'patch.m'; some useful options include 'EdgeColor', 'EdgeAlpha',
%       'LineWidth', and so on. Besides, you can determine whether to use
%       the map projection or not when plotting ('projection', 'on').
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 14 Dec 2024.
% Email: wwu@vims.edu
%
% See also: tripcolor and patch

%% Parse inputs
if numel(find(size(varData)~=1))~=1
    error('the input variable must be a vector!')
end
if numel(varData)~=Mobj.nNodes && numel(varData)~=Mobj.nElems && numel(varData)~=Mobj.nEdges
    error('the length of input variable is not consistent with the given mesh grid!')
end
varData = double(varData(:));

ind = find(strncmpi(varargin, 'projection', 4), 1, 'last' ); % only the last one is valid
if isempty(ind)
    proj_flag = 'off';
else
    proj_flag = varargin{ind+1};
    varargin(ind:ind+1) = [];
end
%% Display
switch proj_flag
    case 'off'    % without projection
        tripcolor(Mobj.tri, Mobj.lon, Mobj.lat, varData, varargin{:})

    case 'on'   % with projection
        proj = projcrs(3395, 'Authority', 'EPSG');  % % EPSG 3395: Mercator projection
        [x, y] = projfwd(proj, Mobj.lat, Mobj.lon);

        tripcolor(Mobj.tri, x, y, varData, varargin{:})
        xticks = get(gca, 'XTick'); yticks = get(gca, 'YTick');

        [~, lon_xticks] = projinv(proj, xticks, mean(y)*ones(size(xticks)));
        [lat_yticks, ~] = projinv(proj, mean(x)*ones(size(yticks)), yticks);

        if fix((max(lon_xticks)-min(lon_xticks))/0.5)>4
            lon_xticks = round(lon_xticks*2)/2;
        else
            lon_xticks = round(lon_xticks,2);
        end

        if fix((max(lat_yticks)-min(lat_yticks))/0.5)>4
            lat_yticks = round(lat_yticks*2)/2;
        else
            lat_yticks = round(lat_yticks,2);
        end

        set(gca, 'XTick', xticks, 'XTickLabel', lon_xticks);
        set(gca, 'YTick', yticks, 'YTickLabel', lat_yticks);
        axis image
end
xlabel('Longitude', 'FontWeight','bold')
ylabel('Latitude', 'FontWeight','bold')
box on
colorbar
colormap(jet(25))
end
