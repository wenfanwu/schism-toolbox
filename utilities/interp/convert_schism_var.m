function varNew = convert_schism_var(Mobj, varRaw, conv_str)
% Convert variables between different mesh centers
%
%% Syntax
% varNew = convert_schism_var(Mobj, varRaw, conv_str)
%
%% Description
% varNew = convert_schism_var(Mobj, varRaw, conv_str)
%
%% Examples 
% varRaw = Mobj.depth;
% varNew = convert_schism_var(Mobj, varRaw, 'node2elem');
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       this datastruct contains mesh info.
% varRaw - raw variable data; double
%       raw variable data @node/element/edge.
% conv_str - conversion string; char
%       conv_str determines the interpolation direction, six different
%       methods are provided now:
%       1) 'node2elem': interpolate from nodes to elements
%       2) 'elem2node': interpolate from elements to nodes
%       3) 'node2edge': interpolate from nodes to edges/sides
%       4) 'edge2node': interpolate from edges to nodes
%       5) 'elem2edge': interpolate from elements to edges (not work yet)
%       6) 'edge2elem': interpolate from edges to elements (not work yet)
%
%% Output Arguments
% varNew - new variable data; double
%       new variable data on required mesh centers.
%
%% Notes
% The function was generated with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 11 Nov 2024. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
switch lower(conv_str)
    case 'node2elem'  % Arithmetic mean at the vertices
        varNew = nan(Mobj.nElems, 1);
        varNew(Mobj.i34==3) = mean(varRaw(Mobj.tri(Mobj.i34==3, 1:3)), 2, 'omitnan');
        varNew(Mobj.i34==4) = mean(varRaw(Mobj.tri(Mobj.i34==4, 1:4)), 2, 'omitnan');

    case 'elem2node'  % Area-weighted average on the surrounding elements
        varTmp = repmat(varRaw(:), [1 4]);
        area = calc_schism_area(Mobj);
        A = repmat(area(:), [1 4]);
        tri = Mobj.tri; tri(isnan(tri)) = Mobj.nNodes+1; % temporally add a node

        varNew = accumarray(tri(:), varTmp(:).*A(:), [], @sum);
        total_weight = accumarray(tri(:), A(:), [], @sum);
        varNew = varNew./total_weight; varNew(end) = [];
        
    case {'node2side', 'node2edge'}  % Arithmetic mean at two endpoints
        varNew = mean(varRaw(Mobj.edg), 2, 'omitnan');

    case {'side2node', 'edge2node'} % Length-weighted average for the surrounding sides.
        varTmp = repmat(varRaw(:), [1 2]);
        [~, edge_len_list, ~] = calc_schism_edge(Mobj);
        L = repmat(edge_len_list(:), [1 2]);
        varNew = accumarray(Mobj.edg(:), varTmp(:).*L(:), [], @sum);

        total_weight = accumarray(Mobj.edg(:), L(:), [], @sum);
        varNew = varNew./total_weight;

    case {'elem2side', 'elem2edge'}
        error('not work yet!')
        
    case {'side2elem', 'edge2elem'}
        error('not work yet!')
end

end