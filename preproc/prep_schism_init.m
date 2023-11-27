function DS = prep_schism_init(Mobj, dst)
% Prepare initial data for SCHISM
%
%% Syntax
% DS = prep_schism_init(Mobj, dst)
%
%% Description
% DS = prep_schism_init(Mobj, dst) prepares the initial data for SCHISM
%
%% Example
% DS = prep_schism_init(Mobj, 'hycom')
%
%% Input Arguments
% Mobj --- the mesh object
% dst --- the data set name;
%
%% Output Arguments
% DS --- a datastruct contains the variable data.
%
%% Notes
% All the variables are stored into the datastruct ('DS') as sub-datastructs,
% thus each variable can have its own independent coordinate system. This
% means that the variables can come from different products.
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
% Last Updated on 2022-05-17.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also: switch

%% Parse inputs
clear DS
switch dst
    case 'hycom_bys'
        DS = get_hycom_init_bys(Mobj);
    case 'hycom_clim'
        DS = get_hycom_init_clim(Mobj);
    case 'hycom_online' % directly download the hycom data from the internet
        varList = {'ssh','temp','salt'};
        hycom_data = get_hycom_online(Mobj.aimpath, Mobj.region, Mobj.time(1), varList);
        for iVar = 1:numel(varList)
            clear D
            varName = varList{iVar};
            D.lon = hycom_data.lon;
            D.lat = hycom_data.lat;
            D.depth = hycom_data.depth;
            D.time = hycom_data.time;
            D.var = hycom_data.(varName);
            DS.(varName) = D;
        end
end
end
