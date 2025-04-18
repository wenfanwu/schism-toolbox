function [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox, margin)
% Extract the index and coordinates of sub-region.
%
%% Syntax
% [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox)
% [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox, margin)
%
%% Description
% [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox) extracts index and
%       coordinates of subregion defined by bbox
% [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox, margin) expands the
%       subregion slightly.
%
%% Examples 
% x = 100:0.25:140;
% y = 20:0.25:45;
% bbox = [110.1, 124.9, 24.8, 30.1];
% [ind_x, ind_y, sx, sy] = sub_region(x, y, bbox);
%
%% Input Arguments
% x - x-axis coordinates; numeric
%       the x-axis coordinates (vector).
% y - y-axis coordinates; numeric
%       the y-axis coordinates (vector).
% bbox - bounding box; numeric
%       the bounding box used to define sub-region. bbox = [x_min, x_max, y_min, y_max];
% margin - the marginal value; numeric
%       the margin value (>0) used to expand the sub-region.
%
%% Output Arguments
% ind_x - index along x-axis; numeric
%       the sub-region index along x-axis.
% ind_y -  index along y-axis; numeric
%       the index along y-axis.
% sx - x-axis coordinates; numeric
%       the x-axis coordinates of sub-region.
% sy - y-axis coordinates; numeric
%       the y-axis coordinates of sub-region.
%
%% Notes
% This function can make sure the sub-region defined by bbox is fully
% covered by the "region" defined by "sx" and "sy". This is important for
% data extra-interpolation near the boundary.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 16 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if nargin < 4; margin = 0; end

min_lon = bbox(1); max_lon = bbox(2);
min_lat = bbox(3); max_lat = bbox(4);

if ~issorted(x, 'ascend') || ~issorted(y, 'ascend')
    error('x and y must be in an ascending order')
end
%% Return the index
margin = max(margin, 0) + 1; % ensure the bbox can be fully covered.
i1 = max(1, find(x >= min_lon, 1, 'first') - margin);
i2 = min(length(x), find(x <= max_lon, 1, 'last') + margin);

j1 = max(1, find(y >= min_lat, 1, 'first') - margin);
j2 = min(length(y), find(y <= max_lat, 1, 'last') + margin);

% index for subregion
ind_x = [i1, i2]; ind_y = [j1, j2];

% subregion x and y
sx = x(min(ind_x):max(ind_x));
sy = y(min(ind_y):max(ind_y));

end