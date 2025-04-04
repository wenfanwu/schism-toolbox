function BdryCnd = interp_schism_bdry(Mobj, DS, varList)
% Interpolate boundary data onto SCHISM vertical layers.
%
%% Syntax
% BdryCnd = interp_schism_bdry(Mobj, DS, varList)
%
%% Description
% BdryCnd = interp_schism_bdry(Mobj, DS, varList) interpolates boundary
% data onto the vertical layers of SCHISM grid.
%
%% Input Arguments
%   Mobj    - Mesh object containing SCHISM mesh and metadata
%   DS      - Struct containing boundary input data (fields: var, depth, time, etc.)
%   varList - Cell array of variable names to interpolate (optional)
%
%% Output Arguments
%   BdryCnd - Struct containing interpolated variables on SCHISM vertical layers
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022. 
% Last Updated on 1 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if nargin < 3
    varList = fieldnames(DS);
end
nVars = numel(varList);

%% Begin to interpolate
for iVar = 1:nVars
    varName = varList{iVar};
    D = DS.(varName);
    
    % Check the index of used open boundary segments 
    % It must be consistent with hgrid files
     [obc_bnds, idx] = ismember(Mobj.obc_nodes(1,:), D.ind); 
     obc_bnds = find(obc_bnds);
    if ~issorted(idx(obc_bnds))
        error('the index of open boundary segments is not consistent with hgrid files')
    end

    % Prepare time index
    model_time = dateshift(Mobj.time, 'start', D.unit_time);  % expand the time
    [idx, ind_full] = ismember(model_time, D.time);
    if sum(idx)~=numel(model_time)
        error('the time range of bdry inputs should cover the model time!')
    end

    % Interpolate the data vertically
    depRaw = abs(D.depth);
    depRaw = depRaw-min(depRaw);  % avoid nan values at 0-m layer
    depBnd = abs(Mobj.depLayers(:, D.ind));

    nt = numel(D.time);
    nps = numel(D.ind);
    switch lower(varName)
        case 'ssh'
            varBnd = D.var;

        otherwise
            varBnd = zeros(nps, Mobj.maxLev, nt);
            for iNode = 1:nps
                varRaw = squeeze(D.var(iNode,:,:));
                varRaw = fillmissing(varRaw, 'previous', 1);  % fill the potential missing values at deep layers
                depNew = fillmissing(depBnd(:, iNode), 'previous');
                varBnd(iNode,:,:) = multi_interp1(depRaw, varRaw, depNew, 1);
            end
    end

    % Kill nan values on each open boundaries
    obc_inds = [0, cumsum(Mobj.obc_lens(obc_bnds))];
    for ii = 1:numel(obc_bnds)
        ind_seg = obc_inds(ii)+1:obc_inds(ii+1);
        var_seg = varBnd(ind_seg,:,:); % data on each open boundary segment
        varBnd(ind_seg,:,:) = kill_nans(var_seg);
    end
    if numel(find(isnan(varBnd(:)))) ~= 0
        error('there are all NaN values in the matrix')
    end

    BdryCnd.(varName) = varBnd(:,:, ind_full);  % expand on the time dimension
end

end

function var_out = kill_nans(var_in)
% Fills NaN values in a 3D matrix using nearest-neighbor logic

if any(isnan(var_in(:)))
    warning on
    warning('NaN values found. Attempting to fill...');
    var_out = fillmissing(var_in, 'nearest', 1);
    var_out = fillmissing(var_out, 'nearest', 2);
    var_out = fillmissing(var_out, 'nearest', 3);
    var_out = fillmissing(var_out, 'constant', mean(var_out(:), 'omitnan'));
else
    var_out = var_in;
end

end















