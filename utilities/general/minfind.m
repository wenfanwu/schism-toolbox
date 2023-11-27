function indMin = minfind(varBase, varFind, topNum)
% Find the position of the closest value
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
% Last Updated on 2022-10-04.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also:

%% Parse inputs
if nargin < 3
    topNum = 1;
end

varBase = varBase(:)';
varFind = varFind(:);

diff_vals = abs(varBase-varFind);
[~, indMin] = sort(diff_vals, 2);
indMin = indMin(:, 1:topNum);

end
