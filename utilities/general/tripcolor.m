function tripcolor(tri, x, y, z, varargin)
% Visualize the variable on an unstructured grid.
%
%% Syntax
% tripcolor(tri, x, y, z, varargin)
% tripcolor(tri, x, y, z)
%
%% Description
% tripcolor(tri, x, y, z) visualizes the variables defined at the
%       centers of nodes/elems/edges.
% tripcolor(tri, x, y, z, varargin) specifies your own options.
%
%% Input Arguments
% tri - connectivity table; double
%       the connectivity table (N*3 or N*4) of the unstructured grid.
% x -  x-axis coordinates; double
%       x-axis coordinates (nNodes*1) defined at node centers.
% y -  y-axis coordinates; double
%       x-axis coordinates  (nNodes*1) defined at node centers.
% z -  variable data; double
%       variable vector (nNodes/nElems/nEdges*1) defined at the centers of
%       nodes/elems/edges.
% varagin - options; char
%       the same as the 'patch' function.
%
%% Notes
% 1) work for mixed triangular/quad grid
% 2) variable z can be defined at the centers of nodes/elems/edges.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 14 Dec 2024.
% Email: wwu@vims.edu
%
% See also: patch

%% Parse inputs
i34 = sum(~isnan(tri), 2);

if numel(find(i34==4)) == 0
    tri = tri(:, 1:3);  % for speeding
end

%% Display
switch length(z(:))
    case numel(x)  % @nodes
        opts = {'FaceColor','interp', 'EdgeAlpha', 0.1, 'CDataMapping','scaled', 'EdgeColor', 'none', 'LineWidth', 0.0025};
        F = fillmissing(tri,"previous",2); % for speeding
        V = [x(:) y(:)]; C = z;

    case size(tri,1)  % @elems
        opts = {'FaceColor','flat', 'EdgeAlpha', 0.1, 'CDataMapping','scaled', 'EdgeColor', 'none', 'LineWidth', 0.0025};
        F = fillmissing(tri,"previous",2); % for speeding
        V = [x(:) y(:)]; C = z;

    otherwise  % @edges
        edg = get_edges(tri);
        nEdges = size(edg,1);

        edg = [edg nan(nEdges,1)]'; % use NaN to sperate all edges
        ind_valid = ~isnan(edg(:));
        edg(isnan(edg)) = 1;
        xe = x(edg(:));
        ye = y(edg(:));
        xe(~ind_valid) = nan;
        ye(~ind_valid) = nan;

        V = [xe(:) ye(:)];
        F = 1:size(V,1);
        C = repmat(z, [1 3])';
        C = C(:); C(~ind_valid) = nan;

        opts = {'EdgeColor','flat','FaceColor','none','LineWidth', 1};
end

varargin = [opts(:)', varargin(:)'];
patch('Faces',F, 'Vertices',V, 'FaceVertexCData',C(:), varargin{:});

% Enable more powerful datatips
dcm = datacursormode; 
dcm.UpdateFcn = @schism_datatips;
end

function edges = get_edges(tri)
% get edge info

nElems = size(tri,1);
if size(tri,2) == 3
    tri = [tri nan(nElems,1)];
end

i34 = sum(~isnan(tri), 2);
ind_e3 = ones(size(tri'));
ind_e3(3, i34==3) = nan;
ind_e3 = ind_e3(:);

tri(i34==3,4) = tri(i34==3,1);

edges = nan(4*nElems, 2);
edges(1:4:end, :) = tri(:, [2, 3]); % start from the 2nd edge
edges(2:4:end, :) = tri(:, [3, 4]);
edges(3:4:end, :) = tri(:, [4, 1]);  % invalid edge for triangular cells
edges(4:4:end, :) = tri(:, [1, 2]);

edges(isnan(ind_e3),:) = [];

% find all the independent edges
sorted_edges = sort(edges, 2);
edges = unique(sorted_edges, 'rows','stable');

end



