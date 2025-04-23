function InitCnd = interp_schism_init(Mobj, DS, varList)
% Interpolate initial data onto SCHISM grids
%
%% Syntax
% InitCnd = interp_schism_init(Mobj, DS)
% InitCnd = interp_schism_init(Mobj, DS, varList)
%
%% Description
% InitCnd = interp_schism_init(Mobj, DS) interpolates initial variables onto SCHISM grids
% InitCnd = interp_schism_init(Mobj, DS, varList) specifies the variables
%
%% Examples 
% DS = prep_schism_init(Mobj, 'hycom_bys'); 
% varList = {'ssh','temp', 'salt'}; 
% InitCnd = interp_schism_init(Mobj, DS, varList);
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the data struct containing mesh info.
% DS - the data struct; datastruct
%       the datastruct containing gridded variables, resulting from "prep_schism_init".
% varList - variable list (optional); cell
%       the variables to be interpolated; e.g., varList = {'ssh','temp', 'salt'};  
%       all available variables in 'DS' will be used by default.
%
%% Output Arguments
% InitCnd - initial data; datastruct
%       the initial data interpolated onto unstructured grids.
%
%% Notes
% 1) this fucntion will fill the NaN values adjacent to the coast based on "inpaint_nans.m".
% 2) this function will fill the NaN values at deeper layers using the bottom-most value
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022. 
% Last Updated on 22 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: interp_schism_bdry

%% Parse inputs
if ~isfield(Mobj, 'depLayers'); error('Please add vertical grids first!'); end
if nargin < 3; varList = {DS.Variable}; end

%% Interpolation
nVars = numel(varList);
InitCnd(nVars, 1) = struct('Variable', [], 'Data', [], 'Time', []);

for iVar = 1:nVars
    varName = varList{iVar};
    disp(['begin to interp the ', varName])
    
    ind_var = strcmp({DS.Variable}, varName);
    D = DS(ind_var);
    lonRaw = D.Lon; latRaw = D.Lat; varRaw = squeeze(D.Data);
    switch dimnum(varRaw)
        case 2
            varNew = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw);
        case 3
            depRaw = abs(D.Depth);  % this will be fixed in the future to consider "negative depth"
            varTmp = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw, [1 1]);
            varNew = interp_deps(depRaw, varTmp, abs(Mobj.depLayers));
        otherwise
            error(['dimension error for ', varName, '!'])
    end
    InitCnd(iVar).Variable = varName;
    InitCnd(iVar).Data = varNew;
    InitCnd(iVar).Time = Mobj.time(1);
end
%% T/S Contraints in SCHISM
ind_var = find(strcmp({InitCnd.Variable}, 'temp'));
if ~isempty(ind_var)
    disp('temp is clipped to [-2, 40] ')
    temp_min = -2; temp_max = 40;
    InitCnd(ind_var).Data = min(temp_max, max(temp_min, InitCnd(ind_var).Data));
end

ind_var = strcmp({InitCnd.Variable}, 'salt');
if ~isempty(ind_var)
    disp('salt is clipped to [0, 42] ')
    salt_min = 0; salt_max = 42;
    InitCnd(ind_var).Data = min(salt_max, max(salt_min, InitCnd(ind_var).Data));
end
end

























































