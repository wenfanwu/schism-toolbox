function varDeps = interp_deps(Mobj, varRec, depRec)
% 
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
%
%
%% Output Arguments
% 
% 
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-26.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if size(varRec, 2) ~= length(depRec)
    varRec = varRec';
end
depRec = abs(depRec);
depRec = depRec-min(depRec);

depLayers = abs(Mobj.depLayers);
varDeps = nan(Mobj.nNodes, Mobj.maxLev);
for iNode = 1:Mobj.nNodes
    progressbar(iNode/Mobj.nNodes)
    
    depTri = depLayers(:, iNode);
    varTmp = varRec(iNode,:);
    varDeps(iNode,:) = multi_interp1(depRec, varTmp, depTri, 1);
end
varDeps = fillmissing(varDeps, 'previous', 2, 'EndValues', 'previous');
end






