function dimNums = dimnum(varData)
% Calculate the # of valid dimensions
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
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-05-18.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if numel(varData) == 1
    dimNums = 0;
else
    dimNums = numel(find(size(varData)~=1));
end

end