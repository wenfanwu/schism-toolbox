function check_schism_init(Mobj, DS, InitCnd, varName)
% Check the interpolation of initial fields (surface/bottom)
%
%% Syntax
% check_schism_init(Mobj, DS, InitCnd, varName)
%
%% Description
% check_schism_init(Mobj, DS, InitCnd, varName) checks the interpolation of initial fields
%
%% Examples 
% check_schism_init(Mobj, DS, InitCnd, 'temp')
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the data struct containing mesh info.
% DS - the data struct; datastruct
%       the datastruct containing gridded variables, resulting from "prep_schism_init".
% InitCnd - initial data; datastruct
%       the initial data on unstructured grids, resulting from "interp_schism_init"
% varName - the variable name; char
%       the variable that needs to be check. e.g., varName = 'temp'
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 1 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: interp_schism_init

%%  Parse inputs
ind_var = strcmp({DS.Variable}, varName); D1 = DS(ind_var);
ind_var = strcmp({InitCnd.Variable}, varName); D2 = InitCnd(ind_var);
switch varName
    case 'ssh'
        varRaw_b = D1.Data; varRaw_s = D1.Data;
        varNew_s = D2.Data; varNew_b = D2.Data;
    otherwise
        varRaw = D1.Data;  varNew = D2.Data;

        [nx, ny, nz] = size(varRaw);
        var2d = reshape(varRaw, nx*ny, nz);
        idx = sum(~isnan(var2d), 2); idx(idx == 0) = 1;  
        linIdx = sub2ind([nx*ny, nz], (1:nx*ny)', idx);

        varRaw_b = reshape(var2d(linIdx), nx, ny);
        varRaw_s = squeeze(varRaw(:,:, 1));

        varNew_b = varNew(sub2ind(size(varNew), sum(~isnan(Mobj.depLayers)), 1:size(varNew,2)));
        varNew_s = varNew(1,:);
end
%% Display
cmap = jet(25);

figure('Color', 'w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile
disp_schism_var(Mobj, varNew_s)
hold on
plot_schism_bnds(Mobj)
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
cb = colorbar;
colormap(cmap)
vm = caxis; %#ok<*CAXIS>
cb.Ticks = linspace(vm(1), vm(2), 6);
cb.TickLabels = compose('%.2f', cb.Ticks);
title(['SCHISM (surface ', varName,')'], 'FontWeight','bold')

nexttile
pcolor(D1.Lon, D1.Lat, varRaw_s')
shading flat
hold on
plot_schism_bnds(Mobj)
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
cb = colorbar;
colormap(cmap);
caxis(vm)
cb.Ticks = linspace(vm(1), vm(2), 6);
cb.TickLabels = compose('%.2f', cb.Ticks);
xlabel('Longitude', 'FontWeight','bold')
ylabel('Latitude', 'FontWeight','bold')
title(['Raw data (surface ', varName,')'], 'FontWeight','bold')

nexttile
disp_schism_var(Mobj, varNew_b)
hold on
plot_schism_bnds(Mobj)
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
cb = colorbar;
colormap(cmap)
vm = caxis; %#ok<*CAXIS>
cb.Ticks = linspace(vm(1), vm(2), 6);
cb.TickLabels = compose('%.2f', cb.Ticks);
title(['SCHISM (bottom ', varName,')'], 'FontWeight','bold')

nexttile
pcolor(D1.Lon, D1.Lat, varRaw_b')
shading flat
hold on
plot_schism_bnds(Mobj)
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
cb = colorbar;
colormap(cmap)
caxis(vm)
cb.Ticks = linspace(vm(1), vm(2), 6);
cb.TickLabels = compose('%.2f', cb.Ticks);
xlabel('Longitude', 'FontWeight','bold')
ylabel('Latitude', 'FontWeight','bold')
title(['Raw data (bottom ', varName,')'], 'FontWeight','bold')

% auto_center
end




