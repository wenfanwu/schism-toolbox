function [Fx, Fy, Fxy] = calc_schism_grad(Mobj, F, grd_ctr)
% Calculate the gradient (units m-1) on an unstructured grid.
%
%% Syntax
% [Fx, Fy, Fxy] = calc_schism_grad(Mobj, F)
%
%% Description
% [Fx, Fy, Fxy] = calc_schism_grad(Mobj, F)
%
%% Examples
% [Fx, Fy, Fxy] = calc_schism_grad(Mobj, Mobj.depth);
%
% figure
% disp_schism_var(Mobj, Fxy, 'EdgeColor', 'k')
% clim([0 0.025])
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
% F - variable data; double
%       variable data defined at node centers.
% grd_ctr - grid center; char
%       the grid center used for output (node/edge). Note that the original
%       calculation is done on element centers. Default: grd_ctr = 'node'.
%
%% Output Arguments
% Fx - gradient along x-axis; double
%       variable gradient along x-axis
% Fy - gradient along y-axis; double
%       variable gradient along y-axis
% Fxy - sum of gradients; double
%       total gradient.
%
%% Notes
% This function calculates the gradient at element centers using Finite
% Volume Method. As for quad cells, the gradient is defined as the average
% of its two splited triangular cells' gradient. This function also
% considered the map projetion.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 8 Jun 2025.
% Email: wwu@vims.edu
%
% See also: calc_schism_area

%% Parse inputs
if strncmpi(Mobj.coord, 'geographic', 3)
    ctype = 1; % coordinate type
else
    ctype = 0;
end

if nargin<3
    grd_ctr = 'node'; 
end

if ctype == 1
    R = 6371000; % earth radius
    lat_rad = (Mobj.lat*2*pi)/360;
    lon_rad = (Mobj.lon*2*pi)/360;

    % convert from lat/lon to distance (m)
    X = R * cos(lat_rad) .* lon_rad;
    Y = R * lat_rad;
else
    X = Mobj.lon;
    Y = Mobj.lat;
end
%% Calculation
[A, ~, A4] = calc_schism_area(Mobj);
i34 = Mobj.i34; F = F(:);

tri = Mobj.tri; tri(Mobj.i34==3, 4) = tri(Mobj.i34==3, 2);
v1 = F(tri(:,1)); v2 = F(tri(:,2)); v3 = F(tri(:,3)); v4 = F(tri(:,4));
x1 = X(tri(:,1)); x2 = X(tri(:,2)); x3 = X(tri(:,3)); x4 = X(tri(:,4));
y1 = Y(tri(:,1)); y2 = Y(tri(:,2)); y3 = Y(tri(:,3)); y4 = Y(tri(:,4));

Fx = nan(Mobj.nElems,1); Fy = nan(Mobj.nElems,1);

% for triangular cells
sx1 = v1.*(y2-y3)+v2.*(y3-y1)+v3.*(y1-y2);
sy1 = (x3-x2).*v1+(x1-x3).*v2+(x2-x1).*v3;
Fx(i34==3) = sx1(i34==3)./(2*A(i34==3));
Fy(i34==3) = sy1(i34==3)./(2*A(i34==3));

% for quad cells
sx2 = v1.*(y3-y4)+v3.*(y4-y1)+v4.*(y1-y3);
sy2 = (x4-x3).*v1+(x1-x4).*v3+(x3-x1).*v4;

% calculation on elements (most accurate)
Fx(i34==4) = (sx1(i34==4)./(2*A4(:,1)) + sx2(i34==4)./(2*A4(:,2)))/2;
Fy(i34==4) = (sy1(i34==4)./(2*A4(:,1)) + sy2(i34==4)./(2*A4(:,2)))/2;
Fxy = hypot(Fx, Fy);

% convert onto node/side centers
switch lower(grd_ctr(1:4))
    case 'node'
        Fx = convert_schism_var(Mobj, Fx, 'elem2node');
        Fy = convert_schism_var(Mobj, Fy, 'elem2node');
        Fxy = convert_schism_var(Mobj, Fxy, 'elem2node');
    case {'edge', 'side'}  % not work yet
        Fx = convert_schism_var(Mobj, Fx, 'elem2edge');
        Fy = convert_schism_var(Mobj, Fy, 'elem2edge');
        Fxy = convert_schism_var(Mobj, Fxy, 'elem2edge');
end

end