function check_schism_icbc(Mobj, varName, ind_lev)
% Check the consistency between hotstart.nc and *.th.nc files
%
%% Syntax
% check_schism_icbc(Mobj, varName)
% check_schism_icbc(Mobj, varName, ind_lev)
%
%% Description
% check_schism_icbc(Mobj, varName) checks the consistency bween I.C. nad B.C.
% check_schism_icbc(Mobj, varName, ind_lev) specifies the vertical layer
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct containing mesh info.
% varName - the variable name; char
%       the variable to be shown; e.g., varName = 'temp'.
%       available variables: temp, salt, ssh/elev.
% ind_lev - the index of level (optional); numeric
%       the index of vertica layer; Default: ind_lev = 1, which means the
%       surface layer.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 1 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: check_schism_init and check_schism_bdry

%% Parse inputs
if nargin < 3; ind_lev = 1; end  % 1 means the first level
idx_time = 1; datapath = Mobj.aimpath;

%% Load data
init_file = fullfile(datapath, 'hotstart.nc');
switch varName
    case {'elev', 'ssh'}
        bdry_file = fullfile(datapath, 'elev2D.th.nc');
        var_init = squeeze(ncread(init_file, 'eta2'));
        ind_lev = Mobj.maxLev;
    case 'temp'
        bdry_file = fullfile(datapath, 'TEM_3D.th.nc');
        var_init = squeeze(ncread(init_file, 'tr_nd', [1 Mobj.maxLev-ind_lev+1 1], [1 1 Mobj.nNodes]));
    case 'salt'
        bdry_file = fullfile(datapath, 'SAL_3D.th.nc');
        var_init = squeeze(ncread(init_file, 'tr_nd', [2 Mobj.maxLev-ind_lev+1 1], [1 1 Mobj.nNodes]));
    otherwise
        disp('not comple now!')
end
varTmp = ncread(bdry_file, 'time_series');  % not complete
var_bnd = squeeze(varTmp(1, Mobj.maxLev-ind_lev+1, :, idx_time));

% Determine the open boundary segments
cum_lens = [0, cumsum(Mobj.obc_lens)];
obc_bnds = 1:(find(cum_lens==size(varTmp,3))-1);
obc_nodes = Mobj.obc_nodes(:, obc_bnds);
obc_nodes(obc_nodes==0) = [];

lon = Mobj.lon(obc_nodes);
lat = Mobj.lat(obc_nodes);
%% Display
figure('Color', 'w');
disp_schism_var(Mobj, var_init)
vm = caxis; %#ok<*CAXIS>
hold on
scatter(lon, lat, 60, var_bnd, 'filled', 'MarkerEdgeColor', [.5 .5 .5], 'MarkerEdgeAlpha', 0.1)
axis image
colormap(jet(25))
box on
caxis(vm)
title(varName, 'FontWeight', 'bold')
% auto_center

figure('Color', 'w')
tiledlayout(2,1,'TileSpacing','tight')
nexttile
% subplot(211)
plot(var_init(obc_nodes), 'LineWidth',2, 'Color','b', 'Marker','.', 'MarkerSize',15)
hold on
plot(var_bnd, 'LineWidth',1, 'Color','r', 'Marker','.', 'MarkerSize',10)
legend({'IC', 'BC'})
box on; grid on
xlabel('Along open boundary nodes', 'FontWeight', 'bold')
ylabel(varName, 'FontWeight', 'bold')

nexttile  
% subplot(212)
plot(var_init(obc_nodes)-var_bnd, 'LineWidth',1, 'Color','m', 'Marker','.', 'MarkerSize', 10)
box on; grid on
xlabel('Along open boundary nodes', 'FontWeight', 'bold')
ylabel([varName, ' (IC - BC)'], 'FontWeight', 'bold')

end













