function S = calc_schism_contour(Mobj, z, levels)
% Extract contour lines from SCHISM grid.
%
%% Syntax
% S = calc_schism_contour(Mobj, z, levels)
%
%% Description
% S = calc_schism_contour(Mobj, z, levels) extracts contour lines from an
%       unstructured grid based on user-specified levels.
%
%% Examples
% levels = [5, 10, 30, 100];
% S = calc_schism_contour(Mobj, Mobj.depth, levels);
%
% figure
% hold on
% for ii = 1:length(S)
%     plot(S(ii).X, S(ii).Y)
% end
%
% shapewrite(S, 'C:\WorkDisk\test')  % save as shapefiles if necessary
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       A datastruct containing the SCHISM mesh information.
%       Required fields:
%         - tri  : element connectivity (Nx4)
%         - lon  : node x-coordinates (nNodes x 1)
%         - lat  : node y-coordinates (nNodes x 1)
%         - i34  : element type (optional; 3 = triangle, 4 = quadrilateral)
% z - scalar variable; double
%       Scalar variable defined at node centers (e.g., depth, salinity).
% levels - contour levels; vector
%       A vector specifying which scalar values to extract contour lines for.
%
%% Output Arguments
% S - shapefile-compatible struct array; datastruct
%       Each struct contains:
%         - Geometry    : 'Line'
%         - X           : x-coordinates of merged contour segments (NaN-separated)
%         - Y           : y-coordinates of merged contour segments
%         - BoundingBox : [minX maxX; minY maxY]
%         - Level       : corresponding contour level
%
%% Notes
% Each element of the output struct array corresponds to one contour level.
% All line segments for that level are merged into a single polyline using
% NaNs as separators. This format is directly compatible with shapewrite()
% for export, or plotting in MATLAB/GIS.
%
% This function supports both triangular and quadrilateral elements. For
% quadrilateral cells, the gradient is computed using all four edges.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 24 Mar 2025.
% Email: wwu@vims.edu
%
% See also: shapewrite

%% Parse inputs
tri = Mobj.tri;
x = Mobj.lon(:);
y = Mobj.lat(:);
z = z(:);

if isfield(Mobj, 'i34')
    i34 = Mobj.i34(:);
else
    i34 = ~isnan(tri(:,4))+3;
end
%% Extract contours
nLevels = length(levels);
S(nLevels, 1) = struct('Geometry', 'Line', 'BoundingBox', [], 'X', [], 'Y', [], 'Level', []);
for iLevel = 1:nLevels
    level = levels(iLevel);
    X_all = [];
    Y_all = [];

    for iNode = 1:size(tri,1)
        if i34(iNode) == 3
            nodes = tri(iNode,1:3);
            edges = [1 2; 2 3; 3 1];
        elseif i34(iNode) == 4
            nodes = tri(iNode,1:4);
            edges = [1 2; 2 3; 3 4; 4 1];
        end

        vx = x(nodes);
        vy = y(nodes);
        vz = z(nodes);
        pts = [];

        for e = 1:size(edges,1)
            i1 = edges(e,1);
            i2 = edges(e,2);
            z1 = vz(i1);
            z2 = vz(i2);

            if ( (z1 < level && z2 > level) || (z1 > level && z2 < level) )
                t = (level - z1) / (z2 - z1);
                px = vx(i1) + t * (vx(i2) - vx(i1));
                py = vy(i1) + t * (vy(i2) - vy(i1));
                pts(end+1,:) = [px, py]; %#ok<AGROW>
            end
        end

        if size(pts,1) == 2
            X_all = [X_all, pts(1,1), pts(2,1), NaN]; %#ok<AGROW>
            Y_all = [Y_all, pts(1,2), pts(2,2), NaN]; %#ok<AGROW>
        end
    end

    % save as datastruct (shapefile-compatible)
    S(iLevel).Geometry = 'Line';
    S(iLevel).BoundingBox = [min(X_all(~isnan(X_all))), max(X_all(~isnan(X_all))); ...
        min(Y_all(~isnan(Y_all))), max(Y_all(~isnan(Y_all)))];
    S(iLevel).X = X_all;
    S(iLevel).Y = Y_all;
    S(iLevel).Level = level;
end

end
