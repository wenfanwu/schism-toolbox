function c_pts = interp_tri(x_pts, y_pts, x_grd, y_grd, c_grd, nan_flags, method)
% Interpolate gridded data onto scattered points
%
%% Syntax
% c_pts = interp_tri(x_pts, y_pts, x_grd, y_grd, c_grd)
% c_pts = interp_tri(x_pts, y_pts, x_grd, y_grd, c_grd, nan_flags)
%
%% Description
% c_pts = interp_tri(x_pts, y_pts, x_grd, y_grd, c_grd) interpolates
%       gridded data onto scattered points
% c_pts = interp_tri(x_pts, y_pts, x_grd, y_grd, c_grd, nan_flags)
%       specifies the nan flags
%
%% Input Arguments
% x_pts - x-coordinates of points; numeric
% y_pts - y-coordinates of points; numeric
% x_grd - x-coordinates of gridded data; numeric
% y_grd - y-coordinates of gridded data; numeric
% c_grd - the gridded data (nx*ny*nz or nx*ny); numeric
% nan_flags - the nan flags (optional); numeric
%       a two-element vector used to determine whether to fill NaN values
%       or not. Default: nan_flags = [1 1], the first '1' will fill NaNs
%       near the coast while the second one will fill NaNs near the bottom.
% method - interpolate method (optional); char
%       the interpolate function; Default: method = 'linear'.
%
%% Output Arguments
% c_pts -  interpolated scattered data.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 15 Apr 2025.
% Email: wwu@vims.edu
%
% See also: interp_deps

%% Parse inputs
if size(c_grd,1) ~= length(x_grd) || size(c_grd,2) ~= length(y_grd)
    error('the gridded data can not match with the provided coordinates');
end
if nargin < 6; nan_flags = [1 1]; end
if nargin < 7; method = 'linear'; end

%% Fill NaN values horizontally/vertically
% Fill missing values at the coast
if nan_flags(1)==1; c_grd(:,:,1) = inpaint_nans(c_grd(:,:,1)); end 
% Fill missing values at the deep layers
if nan_flags(2)==1; c_grd = fillmissing(c_grd, 'previous', 3, 'EndValues', 'previous'); end

%% Begin to interp
% Pre-allocate the matrix
nps = length(x_pts); nz  = size(c_grd, 3);
c_pts = nan(nps, nz);

% Pre-set the interpolate function
F = griddedInterpolant();  
F.GridVectors = {x_grd, y_grd}; 
F.Method = method;
F.ExtrapolationMethod = 'none';
for iz = 1:nz
    F.Values = c_grd(:,:,iz);
    c_pts(:, iz) = F(x_pts, y_pts);
end

% Avoid extrapolating outliers
is_nan = isnan(c_pts);
min_val = min(c_grd(:), [],'omitnan'); 
max_val = max(c_grd(:),[], 'omitnan');
c_pts = min(max(c_pts, min_val), max_val);
c_pts(is_nan) = nan;

end
