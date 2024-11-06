function [lonReg, lon_flag] = check_lons(lonReg, lonAll)
% Check the consistency of longitude vectors
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
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 9 Oct 2024.
% Email: wwu@vims.edu
%
% See also:

lon_flag = 0;
if max(lonAll)>180 && min(lonReg)<0
    lon_flag = -1;
    lonReg(lonReg<0) = lonReg(lonReg<0)+360;
end
if min(lonAll)<0 && max(lonReg)>180
    lon_flag = 1;
    lonReg(lonReg>180) = lonReg(lonReg>180)-360;
end

end