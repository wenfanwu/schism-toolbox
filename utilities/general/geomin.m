function index = geomin(X, Y, x0, y0, N, ctype)
% Find the indices of closet points
%
%% Syntax
% index = geomin(X, Y, x0, y0)
% index = geomin(X, Y, x0, y0, N)
% index = geomin(X, Y, x0, y0, N, ctype)
%
%% Description
% index = geomin(X, Y, x0, y0) finds the index of the closet points.
% index = geomin(X, Y, x0, y0, N) specifies the # of the closest points.
% index = geomin(X, Y, x0, y0, N, ctype) specifies the coordinate type.
%
%% Examples 
% X = 100:0.1:120; Y = 30:0.1:50;
% x0 = [110, 115.2]; y0 = [33, 38];
% index = geomin(X, Y, x0, y0);
%
%% Input Arguments
% X - x-coordinate; double
%       the searching pool of x-coordinates.
% Y - y-coordinate; double
%       the searching pool of y-coordinates.
% x0 - target x-coordinate; double
%       the target x-coordinates to be searched.
% y0 - target y-coordinate; double
%       the target y-coordinates to be searched.
% N - the # of closest points; double
%       the returned # of closest points. default: N = 1;
% ctype - coordinate type; double/char
%       the coordinate types (geographic/cartesian or 1/0); 
%       default ctype = 'geographic'
%
%% Output Arguments
% index - the index; double
%       the index of closed points;
%
%% Notes
% This function was created with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 20 Feb 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if nargin < 5; N = 1; end
if nargin < 6; ctype = 'geographic'; end
if length(X) ~= length(Y); error('Input vectors must have the same length.'); end

N = min(N, length(X));
%% Calculate
switch lower(ctype)
    case {'cartesian', 0}
        distances = hypot(X(:)-x0(:)', Y(:)-y0(:)');

    case {'geographic', 1}
        % from degree to rad
        lat1 = Y(:)*pi/180;  % N×1
        lon1 = X(:)*pi/180; % N×1
        lat2 = y0(:)'*pi/180;  % 1×M
        lon2 = x0(:)'*pi/180;  % 1×M

        % Haversine formula (vectorized)
        dlat = lat2 - lat1;  % N×M
        dlon = lon2 - lon1;  % N×M
        a = sin(dlat/2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        distances = 6371000 * c; 
end

% sort distances for each target point and get the N nearest indices
[~, sorted_indices] = sort(distances, 1);  % sort along the first dimension (across points)

% return the indices of the N closest points for each target
index = sorted_indices(1:N,:);

end