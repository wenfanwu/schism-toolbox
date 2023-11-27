function varNew = interp_schism_elem2node(Mobj, varRaw, method_flag)
% Interpolate variables from @elems onto @nodes
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
% Last Updated on 2022-06-22.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also:

%% Parse inputs
if nargin <3
    method_flag = 1;
end

switch method_flag
    case 1
        varTmp = repmat(varRaw, [1 3]);
        areas = calc_schism_area(Mobj);
        A = repmat(areas, [1 3]);
        varNew = accumarray(Mobj.tri(:), varTmp(:).*A(:), [], @sum);
        wt_sum = accumarray(Mobj.tri(:), A(:), [], @sum);
        varNew = varNew./wt_sum;
end
end