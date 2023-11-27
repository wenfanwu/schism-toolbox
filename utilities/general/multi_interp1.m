function varNew = multi_interp1(timeRaw, varRaw, timeNew, dimTime, varargin)
%  Excecute interp1 along one dimension of a 2-D matrix
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
% Last Updated on 2022-10-21.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
timeRaw = timeRaw(:);
dim_list = setdiff(1:ndims(varRaw), dimTime);
varRawUsed = permute(varRaw, [dimTime, dim_list]);
varNew = interp1(timeRaw, varRawUsed, timeNew, varargin{:});

indRestore = nan(1,ndims(varNew));
indRestore(dimTime) = 1;
indRestore([dim_list]) = setdiff(1:ndims(varNew), 1);

varNew = permute(varNew, indRestore);

end