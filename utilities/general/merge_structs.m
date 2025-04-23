function S_out = merge_structs(varargin)
% Merges multiple datastructs.
%
%% Syntax
% S_out = merge_structs(S1, S2)
% S_out = merge_structs(S1, S2, S3, ...)
%
%% Description
% merge_structs(S1, S2, ...) merges any number of input structure arrays
% (S1, S2, ...) into a unified structure array. Only fields that appear in
% at least one of the inputs are preserved. The output structure has a
% consistent set of fields across all elements, and the fields are ordered
% by priority:
%
% {'Variable', 'Data', 'Lon', 'Lat', 'Depth', 'Nodes', 'Time'}
%
% Missing fields are filled with empty arrays ([]).
%
%% Examples
% S_out = merge_structs(S1, S2)
% Merges two structure arrays with potentially different field sets into one
% with unified and prioritized field order.
%
%% Input Arguments
% S1, S2, ... - structure arrays; struct
%       Structures that may contain a subset of the following fields:
%       'Variable', 'Data', 'Lon', 'Lat', 'Depth', 'Nodes', 'Time'
%
%% Output Arguments
% S_out - merged structure array; struct
%       Combined structure array with aligned and ordered fields. Fields not
%       present in a structure are filled with [].
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 22 Apr 2025.
% Email: wwu@vims.edu
%
% See also: struct, fieldnames, orderfields

%% Parse inputs
% Define field priority order
priorityFields = {'Variable', 'Data', 'Lon', 'Lat', 'Depth', 'Nodes', 'Time'};

% Combine all input structures into one array
allStructs = vertcat(varargin{:});
allStructs = reshape(allStructs, [], 1);

% Identify all fields that actually appear
usedFields = {};
for i = 1:length(allStructs)
    usedFields = union(usedFields, fieldnames(allStructs(i)), 'stable');
end

% Keep only fields that are both used and in the priority list
finalFields = intersect(priorityFields, usedFields, 'stable');

% Fill missing fields with []
for i = 1:length(allStructs)
    for j = 1:length(finalFields)
        f = finalFields{j};
        if ~isfield(allStructs(i), f)
            allStructs(i).(f) = [];
        end
    end
end

% Sort fields by priority
S_out = orderfields(allStructs, finalFields);

end
