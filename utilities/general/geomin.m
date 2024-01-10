function indMin = geomin(lon, lat, lon_site, lat_site)
% Locate the ordinal number of points
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
% Last Updated on 2022-10-22.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if length(lon_site) > 1 && length(lon_site) == length(lat_site)
    indMin = arrayfun(@(x) geomin(lon, lat, lon_site(x), lat_site(x)), 1:length(lon_site));
else
    for iRange = 0.0001:0.01:0.2
        indNear = abs(lon-lon_site)<iRange & abs(lat-lat_site)<iRange;
        if sum(indNear) ~= 0
            lon(~indNear) = nan;
            lat(~indNear) = nan;
            break
        end
    end
    indMin = minfind(hypot(lon,lat), hypot(lon_site, lat_site));
end
end