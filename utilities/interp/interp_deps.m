function varDeps = interp_deps(depRec, varRec, depLayers)
% Interpolate data from z-layers onto sigma-layers
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
[nNodes, nLevs] = size(varRec);  % # of points; # of z-layers
maxLev = size(depLayers, 1); % # of sigma-layers

%% Check inputs
if nLevs ~= length(depRec)
    if nNodes == length(depRec)
        varRec = varRec';
    else
        error('size mismatch between varRec and depRec.');
    end
end
depRec = abs(depRec) - min(abs(depRec));  % Make sure depth starts from zero
depLayers = abs(depLayers);  % Negative depth is not supported so far

%% Interpolate vertically
varDeps = nan(nNodes, maxLev);
for iNode = 1:nNodes
    progressbar(iNode/nNodes)
    
    depTri = depLayers(:, iNode);
    varTmp = varRec(iNode,:);
    varDeps(iNode,:) = interp1(depRec, varTmp, depTri);
end

% Fill NaN values near the bottom
varDeps = fillmissing(varDeps, 'previous', 2, 'EndValues', 'previous');

end






