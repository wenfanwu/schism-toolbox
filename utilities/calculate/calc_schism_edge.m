function [L_cells, L_edges, tri_edg] = calc_schism_edge(Mobj)
% Calculate the side lengths (m) of each element/cell.
%
%% Syntax
% [L_cells, L_edges, tri_edg] = calc_schism_edge(Mobj)
%
%% Description
% [L_cells, L_edges, tri_edg] = calc_schism_edge(Mobj) calculates side lengths.
%
%% Example
% [L_cells, L_edges, tri_edg] = calc_schism_edge(Mobj)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct used to store the mesh info.
%
%% Output Arguments
% L_cells - edge lengths of each cell; numeric
%       L_cells (nElems*4) stores the edge lengths of each element/cell.
% L_edges - length of each edge; numeric
%       L_edges (nEdges*1) stores the length of each edge/side.
% tri_edg - connectivity table; numeric
%       the connectivity table (nElems*4) between elements and edges.
%
%% Author Info
% Created by Wenfan Wu, Virginia Insitute of Marine Science in 2022.
% Last Updated on 10 Oct 2025.
% Email: wwu@vims.edu
%
% See also: calc_schism_reso

%% Parse inputs
if strncmpi(Mobj.coord, 'geographic', 3)
    ntype = 1;
else
    ntype = 0;
end

%% Calculation
tri3 = Mobj.tri(Mobj.i34==3, 1:3);
tri4 = Mobj.tri(Mobj.i34==4, 1:4);

L_cells = nan(Mobj.nElems,4);
% calculate the triangular areas
x1 = Mobj.lon(tri3(:,[1 2])); y1 = Mobj.lat(tri3(:,[1 2]));
x2 = Mobj.lon(tri3(:,[2 3])); y2 = Mobj.lat(tri3(:,[2 3]));
x3 = Mobj.lon(tri3(:,[3 1])); y3 = Mobj.lat(tri3(:,[3 1]));
e1 = calc_edge_lens(x1, y1, ntype);
e2 = calc_edge_lens(x2, y2, ntype);
e3 = calc_edge_lens(x3, y3, ntype);
L_cells(Mobj.i34==3, 1:3) = [e1(:) e2(:) e3(:)];

% calculate the quadrangular areas (split into two triangles)
x1 = Mobj.lon(tri4(:,[1 2])); y1 = Mobj.lat(tri4(:,[1 2]));
x2 = Mobj.lon(tri4(:,[2 3])); y2 = Mobj.lat(tri4(:,[2 3]));
x3 = Mobj.lon(tri4(:,[3 4])); y3 = Mobj.lat(tri4(:,[3 4]));
x4 = Mobj.lon(tri4(:,[4 1])); y4 = Mobj.lat(tri4(:,[4 1]));
e1 = calc_edge_lens(x1, y1, ntype);
e2 = calc_edge_lens(x2, y2, ntype);
e3 = calc_edge_lens(x3, y3, ntype);
e4 = calc_edge_lens(x4, y4, ntype);
L_cells(Mobj.i34==4, 1:4) = [e1(:) e2(:) e3(:) e4(:)];

% side lengths by order
x0 = Mobj.lon(Mobj.edg); y0 = Mobj.lat(Mobj.edg);
L_edges = calc_edge_lens(x0, y0, ntype);

%% Find the connectivity info between elements and sides
tri_edg = nan(Mobj.nElems, 4);  % connectivity table between element and edge/side
edge_c1 = Mobj.edg(:,1)+Mobj.edg(:,2)*1i;
edge_c2 = Mobj.edg(:,2)+Mobj.edg(:,1)*1i;

for ii = 1:3
    switch ii
        case 1
            et = tri3(:,1)+tri3(:,2)*1i;
        case 2
            et = tri3(:,2)+tri3(:,3)*1i;
        case 3
            et = tri3(:,3)+tri3(:,1)*1i;
    end
    [~, idx1] = ismember(et, edge_c1);
    [~, idx2] = ismember(et, edge_c2);
    tri_edg(Mobj.i34==3,ii) = idx1+idx2;
end

for ii = 1:4
    switch ii
        case 1
            et = tri4(:,1)+tri4(:,2)*1i;
        case 2
            et = tri4(:,2)+tri4(:,3)*1i;
        case 3
            et = tri4(:,3)+tri4(:,4)*1i;
        case 4
            et = tri4(:,4)+tri4(:,1)*1i;
    end
    [~, idx1] = ismember(et, edge_c1);
    [~, idx2] = ismember(et, edge_c2);
    tri_edg(Mobj.i34==4,ii) = idx1+idx2;
end

end

function edge_lens = calc_edge_lens(x, y, ntype)
% calculate the edge lengthA

switch ntype
    case 0   % Euclidean distance
        edge_lens = hypot(diff(x,1,2), diff(y,1,2));

    case 1  % spherical distance
        R = 6371000;  % earth radius (m)
        lon_rad = x * pi / 180;
        lat_rad = y * pi / 180;
        dlon = lon_rad(:,2) - lon_rad(:,1);
        dlat = lat_rad(:,2) - lat_rad(:,1);
        a = sin(dlat/2).^2 + cos(lat_rad(:,1)) .* cos(lat_rad(:,2)) .* sin(dlon/2).^2; %  Haversine formula
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        edge_lens = R * c;
end
end