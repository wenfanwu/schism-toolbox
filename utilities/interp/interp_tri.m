function varTri = interp_tri(lonTri, latTri, lonGrid, latGrid, varGrid)
% Interp the gridded data onto triangular
%
%% Syntax
% varTri = interp_tri(lonTri, latTri, lonRaw, latRaw, varRaw)
%
%% Description
% varTri = interp_tri(lonTri, latTri, lonRaw, latRaw, varRaw) interps the
% gridded data onto triangular mesh.
%
%% Example
% varTri = interp_tri(lonTri, latTri, lonRaw, latRaw, varRaw)
% 
%% Input Arguments
% lonTri --- longitue of  trianular nodes
% latTri --- latitude of  trianular nodes
%
%% Notes
% Refine the raw gridded data set first, and then find the nearest value to
% the trianlge grid points. This method is much faster than using scatterInterpolant. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: interp

%% Parse inputs
varGrid = double(varGrid);  % lon*lat*depth, all dimension vector must be in an ascending order.
nVerts = length(lonTri);
nDeps = size(varGrid, 3);

% Fill the missing values at deep layers with the non-nan values at the
% upper layers, to ensure that there are no NANs in the interpolated
% fields. Since the water properties are generally uniform in the deep
% ocean.
varGrid = fillmissing(varGrid, 'previous', 3, 'EndValues', 'previous');

%% Refine the gridded Data
reso_val = 0.02;
[lonRef, latRef, varRef] = refine_data(lonGrid, latGrid, varGrid, reso_val);  % time-consuming

%% Interp onto the triangular mesh
indLon = minfind(lonRef, lonTri);
indLat = minfind(latRef, latTri);

varTri = zeros(nVerts, nDeps);
for iDep = 1:nDeps
    ind_vals = sub2ind(size(varRef), indLon, indLat, iDep*ones(nVerts, 1));
    varTri(:,iDep) = varRef(ind_vals);
end
end

function [lonRef, latRef, varRef] = refine_data(lonRaw, latRaw, varRaw, reso_val)
% Refine the gridded data

lonRef = (lonRaw(1):reso_val:lonRaw(end))';
latRef = (latRaw(1):reso_val:latRaw(end))';

[X, Y] = ndgrid(lonRaw, latRaw);
[X2, Y2] = ndgrid(lonRef, latRef);

nDeps = size(varRaw,3);
varRef = zeros(length(lonRef), length(latRef), nDeps);
for iDep = 1:nDeps
    V = squeeze(varRaw(:,:,iDep));
    V = inpaint_nans(V);   % use inpaint_nans to expand the value in coastal regions.
    F = griddedInterpolant(X, Y, V);
    varRef(:,:,iDep) = F(X2, Y2);
end
end















