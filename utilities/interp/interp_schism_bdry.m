function BdryCnd = interp_schism_bdry(Mobj, DS, varList, bdry_time)
% Interpolate boundary data onto SCHISM vertical layers.
%
%% Syntax
% BdryCnd = interp_schism_bdry(Mobj, DS)
% BdryCnd = interp_schism_bdry(Mobj, DS, varList)
% BdryCnd = interp_schism_bdry(Mobj, DS, varList, bdry_time)
%
%% Description
% BdryCnd = interp_schism_bdry(Mobj, DS) interpolates boundary data onto
%       the vertical layers of SCHISM grid.
% BdryCnd = interp_schism_bdry(Mobj, DS, varList) specifies the variables
%        to be interpolated. 
% BdryCnd = interp_schism_bdry(Mobj, DS, varList, bdry_time) specifies the
%        time series for intertpolation at open boundaries. 
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the data struct used to store mesh info.
% DS - data struct; datastruct
%       the datastruct containing raw boundary input data (fields: var,
%       depth, time, etc.).
% varList - variable list (optional); cell
%       the variable to be interpolated. Default: varList = fieldnames(DS).
% bdry_time - interpolation time; datetime/cell
%       the time series for intertpolation at open boundaries. bdry_time
%       can be specified as "cell", which means that the time resolution
%       vary from variables. Default: bdry_time = Mobj.time; 
%
%% Output Arguments
% BdryCnd - boundry data; datastruct
%       the datastruct containing interpolated variables on SCHISM vertical
%       layers.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022. 
% Last Updated on 22 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: interp_schism_init

%% Parse inputs
if nargin < 3; varList = {DS.Variable}; end
nVars = numel(varList);
if nargin < 4; bdry_time = repmat({Mobj.time}, [nVars, 1]); end
if isdatetime(bdry_time); bdry_time = repmat({bdry_time}, [nVars, 1]); end

%% Begin to interpolate
BdryCnd(nVars, 1) = struct('Variable', [], 'Data', [], 'Depth', [], 'Nodes', [], 'Time', []);

for iVar = 1:nVars
    varName = varList{iVar};
    disp(['begin to interp the ', varName])

    ind_var = strcmp({DS.Variable}, varName);
    D = DS(ind_var);
    % =============== Chek the inputs ================
    % check the open boundary nodes
    [obc_bnds, idx] = ismember(Mobj.obc_nodes(1,:), D.Nodes);
    obc_bnds = find(obc_bnds);
    is_bnd = all(ismember(D.Nodes(:), Mobj.obc_nodes_tot(:)));  % for open boundary nodes
    if ~issorted(idx(obc_bnds)) && is_bnd
        error('the index of open boundary segments is not consistent with the hgrid file')
    end

    % check the time range
    if max(bdry_time{iVar}) > max(D.Time) || min(bdry_time{iVar}) < min(D.Time)
        error('the time range cannot cover the model time!')
    end

    % =========== Interpolate the data vertically ===========
    depRaw = abs(D.Depth); depRaw = depRaw-min(depRaw);  % avoid nan values at 0-m layer
    depBnd = abs(Mobj.depLayers(:, D.Nodes));
    switch lower(varName)
        case 'ssh'
            varBnd = D.Data;
        otherwise
            varBnd = interp_deps(depRaw, D.Data, depBnd);
    end

    % ========== Kill NaN values along each open boundary ==========
    if is_bnd  % for open boundary nodes
        obc_inds = [0; cumsum(Mobj.obc_lens(obc_bnds(:)))];
        for ii = 1:numel(obc_bnds)
            ind_seg = obc_inds(ii)+1:obc_inds(ii+1);
            var_seg = varBnd(:,ind_seg,:); % data on each open boundary segment
            varBnd(:,ind_seg,:) = kill_nans(var_seg);
        end
    end
    if any(isnan(varBnd(:)))
        warning on; 
        warning(['NaN values were found in the variable "', varName,'" !']); 
    end

    BdryCnd(iVar).Variable = varName;
    BdryCnd(iVar).Data = multi_interp1(D.Time, varBnd, bdry_time{iVar}, min(3, ndims(varBnd)));  % interpolate along the time dimension
    BdryCnd(iVar).Depth = depBnd(1:size(varBnd,1),:);
    BdryCnd(iVar).Nodes = D.Nodes(:);
    BdryCnd(iVar).Time = bdry_time{iVar}(:);
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















