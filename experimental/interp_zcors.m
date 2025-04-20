function varNew = interp_zcors(depLayers, varTri, depNew)
% Interp SCHISM outputs onto standard z levels
% 
%% Syntax
% varNew = interp_zcors(depLayers, varTri, depNew)
%
%% Description 
% varNew = interp_zcors(depLayers, varTri, depNew) interpolates the model
% outputs onto standard z levels given by depNew
%
%% Examples
% varNew = interp_zcors(Mobj.depLayers, varTri, -5)
%
%% Input Arguments
% depLayers --- depth layers (nNodes*nLayers) in negative values
% varTri --- the variable on the model mesh
% depNew --- a specified standard depth. e.g. depNew = -5;
%
%% Output Arguments
% varNew --- interpolated values on standard vertical layer
% 
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: inter_tri

%% Parse inputs
if size(varTri,1)~=size(depLayers,1) || size(varTri,2)~=size(depLayers,2)
    error('check the input matrix!')
end
nNodes = size(varTri,2);

varNew = nan(1, nNodes);
% varTri(abs(varTri)>100) = nan;

for iNode = 1:nNodes
    varRaw = varTri(:,iNode);
    depRaw = depLayers(:,iNode);
    indNaN = find(isnan(depRaw), 1);
    if isempty(indNaN)
        indNaN = numel(depRaw)+1;
    end
    if depNew>min(depRaw)
        varNew(iNode) = interp1(depRaw(1:indNaN-1), varRaw(1:indNaN-1), depNew);
    end
end
end


