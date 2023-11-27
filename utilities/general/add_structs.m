function C = add_structs(A, B)
% Merge two structs
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
M = [fieldnames(A)' fieldnames(B)'; struct2cell(A)' struct2cell(B)'];
C = struct(M{:});

end