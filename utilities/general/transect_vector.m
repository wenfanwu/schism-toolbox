function [tvec, nvec] = transect_vector(x, y, side)
% Compute tangent and normal unit vectors along a transect
%
%% Syntax
% [tvec, nvec] = transect_vector(x, y)
% [tvec, nvec] = transect_vector(x, y, side)
%
%% Description
% [tvec, nvec] = transect_vector(x, y) computes the tangent/normal unit vector
% [tvec, nvec] = transect_vector(x, y, side) defines the positve normal direction (lef/right)
%
%% Example
% x = (0:0.2:10)'; y = sin(x);
% [tvec, nvec] = transect_vector(x, y, 'left');
% 
% figure
% plot(x, y, 'LineWidth', 3, 'Marker', '.', 'Color', 'g')
% hold on
% quiver(x, y, tvec(:,1), tvec(:,2), 0.5, 'r'); hold on
% quiver(x, y, nvec(:,1), nvec(:,2), 0.5, 'b');
% legend('Transect', 'Tangent','Normal'); axis equal
%
%% Input Arguments
% x - x-axis coordinate; numeric
%       the x-axis coordinate (N×1);
% y - y-axis coordinate; numeric
%       the y-axis coordinate (N×1);
% side - positive normal direction; char
%       viewing from the first point towards the last, left/right defines
%       the positive normal direction; Default: side = 'left'. 
% 
%% Output Arguments
% tvec - tangent unit vector; numeric
%       N×2 tangent unit vectors [tx ty]
% nvec - normal unit vector; numeric
%       N×2 normal unit vectors  [nx ny]
%
%% Tips
% ut = u(:) .* tvec(:,1) + v(:) .* tvec(:,2);  % Tangent component
% un = u(:) .* nvec(:,1) + v(:) .* nvec(:,2);  % Normal component
%
%% Notes
% This function was created with the help of ChatGPT.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 23 Sep 2025.
% Email: wwu@vims.edu
%
% See also: def_schism_transect

%% Parse inputs
if nargin < 3; side = 'left'; end

%% Tangent vectors (the same length as x,y)
dx = gradient(x(:)); 
dy = gradient(y(:));
tvec = [dx dy];
tvec = tvec ./ max(hypot(tvec(:,1), tvec(:,2)), eps);   % normalize

%% Normal vectors (rotate tangent by ±90 degrees)
switch side
    case 'left'   % left normal = [-ty, tx]
        nvec = [-tvec(:,2), tvec(:,1)];
    case 'right'  % right normal = [ty, -tx]
        nvec = [ tvec(:,2),-tvec(:,1)];
end

end


