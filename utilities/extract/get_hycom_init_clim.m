function DS = get_hycom_init_clim(Mobj)
% Extract climatological hycom data as initial field
% 
%% Syntax
% DS = get_hycom_init_clim(Mobj)
%
%% Description 
% DS = get_hycom_init_clim(Mobj)
%
%% Examples
% DS = get_hycom_init_clim(Mobj)
%
%% Input Arguments
% Mobj --- the mesh object
%
%% Output Arguments
% DS --- a datastruct containing original variables
% 
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
filepath = 'example_hycom_clim_1995_2020.mat'; 
load(filepath) %#ok<LOAD>

cmon = month(Mobj.time(1));

head_info.lon = hycom_clim.lon;
head_info.lat = hycom_clim.lat;
head_info.depth = abs(hycom_clim.depth);
head_info.time = Mobj.time(1);

DS.ssh = head_info;
DS.ssh.var = squeeze(hycom_clim.ssh(:,:,cmon))';

DS.salt = head_info;
DS.salt.var = permute(squeeze(hycom_clim.salinity(:,:,:,cmon)), [2 1 3]);  % lon*lat*depth

DS.temp = head_info;
DS.temp.var = permute(squeeze(hycom_clim.temperature(:,:,:,cmon)), [2 1 3]);

end