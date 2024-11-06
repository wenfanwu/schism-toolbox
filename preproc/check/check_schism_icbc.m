function check_schism_icbc(Mobj, varName, ind_lev)
% Check the consistency for the hotstart.nc and *.th.nc files
%
%% Syntax
% check_schism_bdry(Mobj, varName)
% check_schism_bdry(Mobj, varName, iTime)
% check_schism_bdry(Mobj, varName, iTime, iDep)
%
%% Description
% check_schism_bdry(Mobj, varName) checks the consistence between the
% SCHISM initial and boundary conditions.
% check_schism_bdry(Mobj, varName, iTime) specifies the timetick
% check_schism_bdry(Mobj, varName, iTime, iDep) specifies the depth layer.
%
%% Example
% check_schism_bdry(Mobj, 'temp')
%
%% Input Arguments
% Mobj --- the mesh object
% varName --- the variable name, e.g. 'temp', 'salt', 'elev'. For more
% options, please refer to the SCHISM NC outputs.
% iTime --- the time step; Default: the first time step
% iDep --- the depth layer; Default: sea surface.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 15 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: importdata and cellfun

%% Parse inputs
if nargin < 2
    varName = 'temp';
end
if nargin < 3
    ind_lev = Mobj.maxLev;  % the last row represents the surface
end
iTime = 1;
datapath = Mobj.aimpath;
%% Load
init_file = fullfile(datapath, 'hotstart.nc');
switch varName
    case 'elev'
        bdry_file = fullfile(datapath, 'elev2D.th.nc');
        varInit = squeeze(ncread(init_file, 'eta2'));
        ind_lev = 1;
    case 'temp'
        bdry_file = fullfile(datapath, 'TEM_3D.th.nc');
        varInit = squeeze(ncread(init_file, 'tr_nd', [1 ind_lev 1], [1 1 Mobj.nNodes]));
    case 'salt'
        bdry_file = fullfile(datapath, 'SAL_3D.th.nc');
        varInit = squeeze(ncread(init_file, 'tr_nd', [2 ind_lev 1], [1 1 Mobj.nNodes]));
    otherwise
        disp('not comple now!')
end
varTmp = ncread(bdry_file, 'time_series');  % not complete
varBnd = squeeze(varTmp(1, ind_lev, :, iTime));

ind = Mobj.obc_nodes_tot;
lon = Mobj.lon(ind);
lat = Mobj.lat(ind);
%% Display
figure('Color', 'w');
disp_schism_var(Mobj, varInit)
varLim = caxis; %#ok<*CAXIS>
hold on
scatter(lon, lat, 60, varBnd, 'filled', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', 0.1)
axis image
colormap(jet(25))
box on
caxis(varLim)
auto_center

end













