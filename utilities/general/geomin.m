function indMin = geomin(lonAll, latAll, siteLon, siteLat, N)
% Find indices of the nearest points on geographic maps
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
% lonAll - longitude vector; double
%       the longitude vector for searching.
% latAll - latitude vector; double
%       the latitude vector for searching.
% siteLon - target longitude; double
%       the target longitude for searching.
% siteLat - target latitude; double
%       the target latitude for searching.
% N - the returned # of points; double
%       the returned # of the first closest points. default: N=1
%
%% Output Arguments
% indMin - index matrix; double
%       the indices of nearest points. the returned indMin is a matrix of
%       (N*M), with M being the length of siteLon/siteLat.  
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 11 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if nargin < 5 
    N = 1;
end
%% Calculate
if max(lonAll)>360
    fcn = @(x,y) hypot(lonAll-x, latAll-y)'; % units
else
    fcn = @(x,y) distance(y,x, latAll, lonAll, [6378.137 0.0818191910428158])';  % km
end

dist_cell = arrayfun(@(x,y) fcn(x,y), siteLon, siteLat, 'UniformOutput',false);
dist = cell2mat(dist_cell);

% sort distances for each target point and get the N nearest indices
[~, sorted_indices] = sort(dist, 2);  % sort along the first dimension (across points)

% return the indices of the N closest points for each target
indMin = sorted_indices(:,1:N); 

end
