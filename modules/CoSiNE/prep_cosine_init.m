function DS = prep_cosine_init(Mobj, dst)
% Prepare initial data for the CoSiNE module
% 
%% Syntax
% DS = prep_cosine_init(Mobj, dst)
%
%% Description 
% DS = prep_cosine_init(Mobj, dst) prepares the initial data for SCHISM
%
%% Example
% DS = prep_cosine_init(Mobj, 'glorys')
%
%% Input Arguments
% Mobj --- the mesh object
% dst --- the data set name
%
%% Output Arguments
% DS --- the datastruct contains the variable data. 
% 
%% Notes
% All the variables are stored into the datastruct ('DS') as sub-datastructs,
% thus each variable can have its own independent coordinate system. This 
% means that the variables can originate from different products.
%
% This is just a simple packaging function, so you can add your own data
% set in this function if needed, just make sure the output format is the
% same.
%
% There are TWO things to note about the format:
% (1) the longitude and latitude vectors must be in an ascending order.
% (2) the dimensions of the variables must be (lon*lat*depth) or (lon*lat)
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2023-11-26.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: switch

%% Parse inputs
switch dst
    case 'test_data'
        load('example_woa_bgc_init_data.mat') %#ok<LOAD>
        DS = woa_init_data;
end

end
