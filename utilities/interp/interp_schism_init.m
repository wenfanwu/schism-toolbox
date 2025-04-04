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
% Last Updated on 1 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: interp_schism_bdry

%% Parse inputs
if ~isfield(Mobj, 'depLayers'); error('Please add vertical grids first!'); end
if nargin < 3; varList = fieldnames(DS); end

%% Interpolation
nVars = numel(varList);
for iVar = 1:nVars
    varName = varList{iVar};
    D = standard_output(DS.(varName)); % ensure the format of D is standard.
    lonRaw = D.lon; latRaw = D.lat; varRaw = squeeze(D.var);

    disp(['begin to interp the ', varName])
    switch dimnum(varRaw)
        case 2
            varNew = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw);
        case 3
            depRaw = abs(D.depth);  % this will be fixed in the future to consider "negative depth"
            varTmp = interp_tri(Mobj.lon, Mobj.lat, lonRaw, latRaw, varRaw);
            varNew = interp_deps(depRaw, varTmp, Mobj.depLayers);
        otherwise
            error(['dimension error for ', varName, '!'])
    end
    InitCnd.(varName) = varNew;
end
%% T/S Contraints in SCHISM
if isfield(InitCnd, 'temp')
    disp('temp is clipped to [-2, 40] ')
    temp_min = -2; temp_max = 40;
    InitCnd.temp = max(temp_min, InitCnd.temp);
    InitCnd.temp = min(temp_max, InitCnd.temp);
end
if isfield(InitCnd, 'salt')
    disp('salt is clipped to [0, 42] ')
    salt_min = 0; salt_max = 42;
    InitCnd.salt = max(salt_min, InitCnd.salt);
    InitCnd.salt = min(salt_max, InitCnd.salt);
end
end

function D = standard_output(D)
% 1) 'lon', 'lat', 'depth' vectors must be in ascending order;
% 2) the region defined by lon/lat must cover you model domain
% 3) 'var' must have the dimensions of lon*lat or lon*lat*depth.

% Ensure "var" is lon x lat (x depth)
if length(D.lon)==length(D.lat)
    warning on
    warning('the lengths of longitude and latitude are the same, be careful')
end
if dimnum(D.var)==2
    if size(D.var,1) == length(D.lat) && size(D.var,2) == length(D.lon)
        D.var = D.var';  % transpose to lon x lat
    end
elseif dimnum(D.var) == 3
    if size(D.var,1) == length(D.lat) && size(D.var,2) == length(D.lon)
        D.var = permute(D.var, [2 1 3]);  % permute to lon x lat x depth
    end
end

% Ensure "lon" is in ascending order
if ~issorted(D.lon); [D.lon, ix] = sort(D.lon); D.var = D.var(ix,:,:); end

% Ensure "lat" is in ascending order
if ~issorted(D.lat); [D.lat, iy] = sort(D.lat); D.var = D.var(:,iy,:); end

% Ensure depth is in ascending order and positive
if dimnum(D.var) == 3 && isfield(D, 'depth') && ~isempty(D.depth)
    if ~issorted(D.depth)
        [D.depth, iz] = sort(D.depth);
        D.var = D.var(:,:,iz);
    end
end
end

























































