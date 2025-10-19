function [S, S3, S4] = calc_schism_area(Mobj)
% Calculate the element area (m^2)
%
%% Syntax
% [S, S3, S4] = calc_schism_area(Mobj)
%
%% Description
% [S, S3, S4] = calc_schism_area(Mobj) calculates the element areas of
%       unstructured grids (units: m^2).  
%
%% Example
% S = calc_schism_area(Mobj);
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
%
%% Output Arguments
% S - element areas; numeric
%       areas of elements  (units: m^2).
% S3 - areas of triangular elements; numeric
%       areas of triangular elements (units: m^2).
% S4 - areas of quadrangular elements; numeric
%       areas of quadrangular elements  (units: m^2).
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 05 Nov 2024. 
% Email: wwu@vims.edu
%
% See also: calc_schism_edge

%% Parse inputs
if strncmpi(Mobj.coord, 'geographic', 3)
    ntype = 1;
else
    ntype = 0;
end

%% Calculation
tri3 = Mobj.tri(Mobj.i34==3, 1:3);
tri4 = Mobj.tri(Mobj.i34==4, 1:4);

% calculate the triangular areas
x1 = Mobj.lon(tri3(:,[1 2])); y1 = Mobj.lat(tri3(:,[1 2]));
x2 = Mobj.lon(tri3(:,[2 3])); y2 = Mobj.lat(tri3(:,[2 3]));
x3 = Mobj.lon(tri3(:,[3 1])); y3 = Mobj.lat(tri3(:,[3 1]));
e1 = calc_edge_lens(x1, y1, ntype);
e2 = calc_edge_lens(x2, y2, ntype);
e3 = calc_edge_lens(x3, y3, ntype);

S3 = calc_triangle_area(e1,e2,e3);

% calculate the quadrangular areas (split into two triangles)
% triangle at (v1,v2,v3)
x1 = Mobj.lon(tri4(:,[1 2])); y1 = Mobj.lat(tri4(:,[1 2]));
x2 = Mobj.lon(tri4(:,[2 3])); y2 = Mobj.lat(tri4(:,[2 3]));
x3 = Mobj.lon(tri4(:,[3 1])); y3 = Mobj.lat(tri4(:,[3 1]));
e1 = calc_edge_lens(x1, y1, ntype);
e2 = calc_edge_lens(x2, y2, ntype);
e3 = calc_edge_lens(x3, y3, ntype);

area4_p1 = calc_triangle_area(e1,e2,e3);

% triangle at (v1,v3,v4)
x1 = Mobj.lon(tri4(:,[1 3])); y1 = Mobj.lat(tri4(:,[1 3]));
x2 = Mobj.lon(tri4(:,[3 4])); y2 = Mobj.lat(tri4(:,[3 4]));
x3 = Mobj.lon(tri4(:,[4 1])); y3 = Mobj.lat(tri4(:,[4 1]));
e1 = calc_edge_lens(x1, y1, ntype);
e2 = calc_edge_lens(x2, y2, ntype);
e3 = calc_edge_lens(x3, y3, ntype);

area4_p2 = calc_triangle_area(e1,e2,e3);

S4 = [area4_p1(:) area4_p2(:)];

S = nan(Mobj.nElems,1);
S(Mobj.i34==3) = S3;
S(Mobj.i34==4) = area4_p1+area4_p2;
end

function area = calc_triangle_area(a,b,c)
% calculate the triangle area using Haron formula

s = (a + b + c) / 2;
area = sqrt(s.*(s-a).*(s-b).*(s-c));
end

function edge_lens = calc_edge_lens(x, y, ntype)
% calculate the edge length using Haversine formula

switch ntype
    case 1
        R = 6371000;  % earth radius (m)
        lon_rad = x * pi / 180;
        lat_rad = y * pi / 180;
        dlon = lon_rad(:,2) - lon_rad(:,1);
        dlat = lat_rad(:,2) - lat_rad(:,1);
        a = sin(dlat/2).^2 + cos(lat_rad(:,1)) .* cos(lat_rad(:,2)) .* sin(dlon/2).^2; 
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        edge_lens = R * c;
    case 0
        edge_lens = hypot(diff(x,1,2), diff(y,1,2));
end
end
