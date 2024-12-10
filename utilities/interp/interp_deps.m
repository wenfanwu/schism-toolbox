function varDeps = interp_deps(depRec, varRec, depLayers)
% interpolate data from standard layers onto sigma layers
% 
%% Syntax
% varDeps = interp_deps(Mobj, varRec, depRec)
%
%% Description 
% varDeps = interp_deps(Mobj, varRec, depRec)
%
%% Input Arguments
% varRec - original variable matrix; double
%       variable matrix (nNodes*nz) on standard layers;
% depRec - standard depths; double
%       standard depth layers;
% depLayers - depth matrix; double
%       the depth matrix (maxLev*nNodes), indicating depth layers for all
%       nodes. It is typically from Mobj.
%
%% Output Arguments
% varDeps - interpolated variable matrix; double
%       varaible matrix (nNodes*maxLev) on sigma layers
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 09 Dec 2024.
% Email: wwu@vims.edu
%
% See also: interp_tri

%% Parse inputs
if size(varRec, 2) ~= length(depRec)
    varRec = varRec';
end
depRec = abs(depRec);
depRec = depRec-min(depRec);

depLayers = abs(depLayers);
[maxLev, nNodes] = size(depLayers);

varDeps = nan(nNodes, maxLev);
for iNode = 1:nNodes
    progressbar(iNode/nNodes)
    
    depTri = depLayers(:, iNode);
    varTmp = varRec(iNode,:);
    varDeps(iNode,:) = interp1(depRec, varTmp, depTri);
end
varDeps = fillmissing(varDeps, 'previous', 2, 'EndValues', 'previous');
end






