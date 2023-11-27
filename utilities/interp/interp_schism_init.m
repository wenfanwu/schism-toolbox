function InitCnd = interp_schism_init(Mobj, DS, varList)
% Interpolate initial data on the mesh grid
% 
%% Syntax
% 
%
%% Description 
% 
%
%% Example
%
%
%% Input Arguments
%
%
%% Output Arguments
% 
% 
%% Tips
% the input option 'varList' enables you to process partial fields in the
% datastruct, this is quite useful when the RAM of your PC is not enough to
% afford a data struct of  too many fields. This case may be encountered
% especially when the model domain is very large.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-05-17.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if ~isfield(Mobj, 'depLayers')
    error('Please add vertical grids first!')
end

if nargin < 3
    varList = fieldnames(DS);
end

%% Interpolation
nVars = numel(varList);
for iVar = 1:nVars
    varName = varList{iVar};
    D = DS.(varName);

    lonRaw = D.lon;
    latRaw = D.lat;
    depRaw = abs(D.depth);
    varRaw = squeeze(D.var);

    disp(['begin to interp the ', varName])
    switch dimnum(varRaw)
        case 2
            varNew = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw);
        case 3
            varTri = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw);
            varNew = interp_deps(Mobj, varTri, depRaw);
        otherwise
            error(['dimension error for ', varName, '!'])
    end
    InitCnd.(varName) = varNew;
end

%% T/S Contraints in SCHISM
if sum(contains(varList, 'temp'))>0
    temp_min = -2; temp_max = 40;
    InitCnd.temp = max(temp_min, InitCnd.temp);
    InitCnd.temp = min(temp_max, InitCnd.temp);
end
if sum(contains(varList, 'salt'))>0
    salt_min = 0; salt_max = 42;
    InitCnd.salt = max(salt_min, InitCnd.salt);
    InitCnd.salt = min(salt_max, InitCnd.salt);
end
end



























































