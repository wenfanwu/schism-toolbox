function var_tri = interp_tri(lon_tri, lat_tri, lon_grid, lat_grid, var_grid)
% Interpolates data from an orthogonal grid onto scattered points
%
%% Syntax
% var_tri = interp_tri_fast(lon_tri, lat_tri, lon_grid, lat_grid, var_grid)
%
%% Description
% var_tri = interp_tri_fast(lon_tri, lat_tri, lon_grid, lat_grid, var_grid)
%       uses the original orthogonal grid to build a gridded interpolant and
%       directly interpolate values at the (lon_tri, lat_tri) locations.
%
%% Input Arguments
% lon_tri, lat_tri: P x 1 scattered x- and y-coordinates
% lon_grid, lat_grid : M x 1 and N x 1 arrays of original grid longitudes and latitudes.
% var_grid : M x N x K array of variable values on the orthogonal grid. K means the depth dimension.
%
%% Output Arguments
%   var_tri : P x K array of interpolated values at each target point and depth level.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 09 Dec 2024.
% Email: wwu@vims.edu
%
% See also: interp_deps

%% Parse inputs
if size(var_grid,1) ~= length(lon_grid) || size(var_grid,2) ~= length(lat_grid)
    error('The first two dimensions of var_grid must match the lengths of lon_grid and lat_grid.');
end

nNodes = length(lon_tri);
nLevs  = size(var_grid, 3);
var_tri = NaN(nNodes, nLevs);
%% Fill NaN values near the bottom
var_grid = fillmissing(var_grid, 'previous', 3, 'EndValues', 'previous');

%% Begin to interp
% If your data contains NaN and you want to fill them, you can use inpaint_nans:
% For each depth level, build the interpolant and perform interpolation

% Create the full grid (M x N)
[LonGrid, LatGrid] = ndgrid(lon_grid, lat_grid);

for iLev = 1:nLevs
    V = var_grid(:,:,iLev);
    V = inpaint_nans(V);  % fill missing values at the coast
    F = griddedInterpolant(LonGrid, LatGrid, V, 'linear', 'none');
    var_tri(:, iLev) = F(lon_tri, lat_tri);
end

% Avoid extrapolating outliers
min_val = min(var_grid(:), [],'omitnan'); 
max_val = max(var_grid(:),[], 'omitnan');
var_tri = min(max(var_tri, min_val), max_val);
end
