function  varNew = get_bnd_vars(lon_bnd, lat_bnd, lon, lat, varRaw)
% Find the variable matrix along the boundary nodes
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
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-10-25.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
nLons = numel(lon);
nLats = numel(lat);

if size(varRaw, 1) ~= nLons
    error('The Longitude dimension is not consistent!')
end
if size(varRaw, 2) ~= nLats
    error('The Latitude dimension is not consistent!')
end

indLon = minfind(lon, lon_bnd);
indLat = minfind(lat, lat_bnd);
indBnd = sub2ind([nLons, nLats], indLon, indLat);

size_list = size(varRaw);
size_tmp = [];
for ii = 3:numel(size_list)
    size_tmp = [size_tmp size(varRaw, ii)]; %#ok<AGROW>
end

varTmp = reshape(varRaw, [size(varRaw,1)*size(varRaw,2) size_tmp]); 
varNew = varTmp(indBnd, :, :);

end