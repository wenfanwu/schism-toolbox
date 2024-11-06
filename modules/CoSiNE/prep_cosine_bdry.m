function DS = prep_cosine_bdry(Mobj, dst)
% Prepare boundary inputs for the CoSiNE module
% 
%% Syntax
% DS = prep_cosine_bdry(Mobj, dst)
%
%% Description 
% DS = prep_cosine_bdry(Mobj, dst) prepares the initial data for SCHISM
%
%% Example
% DS = prep_cosine_bdry(Mobj, 'test_data')
%
%% Input Arguments
% Mobj --- the mesh object
% dst --- the data set name
%
%% Output Arguments
% DS --- the datastruct contains the variable data. 
% 
%% Notes
%
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
        load('example_fake_bgc_data_for_bdry.mat') %#ok<LOAD> 
        DS = fake_bdry_data;
end

end