function varTri = interp_deps(depRaw, varRaw, depTri)
% Interpolate data from z-layers onto sigma-layers
% 
%% Syntax
% varTri = interp_deps(depRaw, varRaw, depTri)
%
%% Description 
% varTri = interp_deps(depRaw, varRaw, depTri) interpolates data from
%       z-layers onto sigma-layers.
%
%% Input Arguments
% depRaw - standard depths; numeric
%       standard depth layers (nz_z*1);
% varRaw - original variable; numeric
%       variable matrix (nz_z*nps*nt or nz_z*nps) on standard layers;
% depTri - depth layers; numeric
%       the depth matrix (nz_s*nps), indicating the depth layers for
%       all nodes.
%
%% Output Arguments
% varTri - interpolated variable; double
%       varaible matrix (nz_s*nps) on sigma layers
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 22 Apr 2025.
% Email: wwu@vims.edu
%
% See also: interp_tri

%% Parse inputs
depRaw = abs(depRaw); depRaw = depRaw-min(depRaw);  % avoid nan values at 0-m layer
depTri = fillmissing(abs(depTri), 'previous', 1, 'EndValues', 'previous');

[nz, nps] = size(depTri);  % nz*nps
nt = 1; if ndims(varRaw)==3; nt = size(varRaw,3); end

varTri = zeros(nz, nps, nt);
for iNode = 1:nps
    varTmp = squeeze(varRaw(:, iNode, :));
    varTmp = fillmissing(varTmp, 'previous', 1, 'EndValues', 'previous');  % fill the potential missing values at deep layers
    depNew = depTri(:, iNode);
    varTri(:,iNode,:) = multi_interp1(depRaw, varTmp, depNew, 1);  % interpolate along the vertical dimension
end

% Fill NaN values near the bottom
varTri = fillmissing(varTri, 'previous', 1, 'EndValues', 'previous');
varTri = squeeze(varTri);

end






