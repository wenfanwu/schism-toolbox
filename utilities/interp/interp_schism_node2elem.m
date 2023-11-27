function varElem = interp_schism_node2elem(Mobj, varNode, method_flag)
% interpolate variables from @node onto @elem
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
if nargin<3
    method_flag = 1;
end

switch method_flag
    case 1
        varElem = mean(varNode(Mobj.tri), 2, 'omitnan');
end

end