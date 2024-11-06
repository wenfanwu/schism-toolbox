function varNew = convert_schism_var(Mobj, varRaw, conv_str)
% Convert variables between different mesh centeres
%
%% Syntax
%
%
%% Description
%
%
%% Examples 
% 
%
%% Input Arguments
% AtmForc - atmospheric forcing; datastruct
%       this datastruct stores atmospheric forcing data.
% suffix_name - filename suffix; char
%       suffix name of the NetCDF files (air/prc/rad).
%
%% Output Arguments
%
%% Notes
% The function was generated with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 24 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: 

switch lower(conv_str)
    case 'node2elem'  % Arithmetic mean at the vertices
        varNew = nan(Mobj.nElems, 1);
        varNew(Mobj.i34==3) = mean(varRaw(Mobj.tri(Mobj.i34==3, 1:3)), 2, 'omitnan');
        varNew(Mobj.i34==4) = mean(varRaw(Mobj.tri(Mobj.i34==4, 1:4)), 2, 'omitnan');

    case 'elem2node'  % Area-weighted average on the surrounding elements
        varTmp = repmat(varRaw(:), [1 4]);
        area = calc_schism_area(Mobj);
        A = repmat(area, [1 4]);
        tri = Mobj.tri;
        tri(isnan(tri)) = Mobj.nNodes+1; % temporally add a node

        varNew = accumarray(tri(:), varTmp(:).*A(:), [], @sum);
        total_weight = accumarray(tri(:), A(:), [], @sum);
        varNew = varNew./total_weight;
        varNew(end) = [];
        
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

    case {'side2elem', 'edge2elem'}
end

end