function BdryCnd = interp_schism_bdry(Mobj, DS, varList)
% Interpolate boundary data onto SCHISM vertical layers.
%
%% Syntax
% BdryCnd = interp_schism_bdry(Mobj, DS, varList)
%
%% Description
% BdryCnd = interp_schism_bdry(Mobj, DS, varList) interpolates boundary
%       data onto the vertical layers of SCHISM grid.
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
% Last Updated on 17 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: interp_schism_init

%% Parse inputs
if nargin < 3; varList = fieldnames(DS); end

%% Begin to interpolate
nVars = numel(varList);
for iVar = 1:nVars
    varName = varList{iVar};
    D = DS.(varName);

    % =============== Chek the inputs ================
    % check the open boundary nodes
     [obc_bnds, idx] = ismember(Mobj.obc_nodes(1,:), D.ind); 
     obc_bnds = find(obc_bnds);
    if ~issorted(idx(obc_bnds))
        error('the index of open boundary segments is not consistent with hgrid files')
    end
    % check the time range
    if max(Mobj.time) > max(D.time) || min(Mobj.time) < min(D.time)
        error('the time range cannot cover the model time!')
    end

    % =========== Interpolate the data vertically ===========
    depRaw = abs(D.depth); depRaw = depRaw-min(depRaw);  % avoid nan values at 0-m layer
    depBnd = abs(Mobj.depLayers(:, D.ind));

    nt = numel(D.time); nps = numel(D.ind);
    switch lower(varName)
        case 'ssh'
            varBnd = D.var;
        otherwise
            varBnd = zeros(nps, Mobj.maxLev, nt);
            for iNode = 1:nps
                varRaw = squeeze(D.var(iNode,:,:));
                varRaw = fillmissing(varRaw, 'previous', 1);  % fill the potential missing values at deep layers
                depNew = fillmissing(depBnd(:, iNode), 'previous');
                varBnd(iNode,:,:) = multi_interp1(depRaw, varRaw, depNew, 1);  % interpolate along the vertical dimension
            end
    end

    % ========== Kill NaN values at each open boundaries ==========
    obc_inds = [0, cumsum(Mobj.obc_lens(obc_bnds))];
    for ii = 1:numel(obc_bnds)
        ind_seg = obc_inds(ii)+1:obc_inds(ii+1);
        var_seg = varBnd(ind_seg,:,:); % data on each open boundary segment
        varBnd(ind_seg,:,:) = kill_nans(var_seg);
    end
    if any(isnan(varBnd(:))); error('NaN values were found!'); end
    
    BdryCnd.(varName) = multi_interp1(D.time, varBnd, Mobj.time, min(3, ndims(varBnd)));  % interpolate along the time dimension
end

end

function var_out = kill_nans(var_in)
% Kill NaNs in a 3-D matrix using nearest-neighbor values

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















