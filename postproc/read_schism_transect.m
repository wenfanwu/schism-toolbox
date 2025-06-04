function [var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri, method)
% Read SCHISM variable along a given transect.
%
%% Syntax
% [var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri)
%
%% Description
% [var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri)
%
%% Examples 
% [var2d, dist2d, dep2d] = read_schism_transect(Mobj, sect_info, var_tri);
% figure; pcolor(dist2d/1e3, dep2d, var2d); shading interp
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% sect_info - section info; datastruct
%       the section datastruct, which is typically created by the
%       function "def_schism_transect". Four fields are required as below:
%       .lon   - Longitude of transect points
%       .lat   - Latitude of transect points
%       .dist  - Cumulative distance along the transect
%       .depth - Depth at each transect point (2D: depth x dist)
% var_tri - variable array; numeric
%       the variable (maxLev*nNodes) on the schism grid.
% method - interpolate method; char
%       the interpolation method; Default: method = 'linear'.
%
%% Output Arguments
% var2d - Interpolated variable along transect (2D: depth x dist)
% dist2d - Distance (meters) along transect (2D: depth x dist)
% dep2d - Depth profile along transect (2D: depth x dist)
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 1 Jun 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
x_pts = sect_info.lon(:);  y_pts = sect_info.lat(:); var_tri = double(var_tri);
assert(numel(x_pts) == numel(y_pts), 'Transect lon/lat must be the same length.');

if nargin < 4; method = 'linear'; end

% Limit interpolation region for speed
buffer = 0.25; 
bbox = [min(x_pts)-buffer, max(x_pts)+buffer, min(y_pts)-buffer, max(y_pts)+buffer];
msk = Mobj.lon > bbox(1) & Mobj.lon < bbox(2) & Mobj.lat > bbox(3) & Mobj.lat < bbox(4);

%% Begin to interpolate
% Interpolate each level separately
nps = numel(x_pts); % # of points
var2d = nan(Mobj.maxLev, nps);
for k = 1:Mobj.maxLev
    data_layer = var_tri(k, msk);
    if all(isnan(data_layer)); continue; end  % skip empty layers
    F = scatteredInterpolant(Mobj.lon(msk), Mobj.lat(msk), data_layer(:), method, 'none');
    var2d(k, :) = F(x_pts, y_pts);
end

% Fill bottom missing values for visualization
dist2d = repmat(sect_info.dist(:)', [Mobj.maxLev, 1]);
dep2d = fillmissing(sect_info.depth, 'previous', 1);
var2d = fillmissing(var2d, 'previous', 1);  % fill vertically

end