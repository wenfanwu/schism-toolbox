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
D = DS.(varName);
switch varName
    case 'ssh'
        varRaw_b = D.var; varRaw_s = D.var;
        varNew_s = InitCnd.ssh; varNew_b = InitCnd.ssh;
    otherwise
        msk_3d = isnan(Mobj.depLayers)';
        InitCnd.(varName)(msk_3d) = nan;

        varRaw = D.var;
        varNew = InitCnd.(varName);

        varRaw_2d = reshape(varRaw, [size(varRaw, 1)*size(varRaw,2) size(varRaw,3)]);
        ind_btm_2d = sum(~isnan(varRaw_2d'));ind_btm_2d(ind_btm_2d==0) = 1;
        varTmp = arrayfun(@(x) varRaw_2d(x, ind_btm_2d(x)), 1:size(varRaw_2d,1));

        varRaw_b = reshape(varTmp, size(varRaw, [1 2]));
        varRaw_s = squeeze(varRaw(:,:, 1));
        
        ind_btm = sum(~isnan(Mobj.depLayers));
        varNew_b = arrayfun(@(x) varNew(x, ind_btm(x)), 1:Mobj.nNodes);
        varNew_s = varNew(:, 1);
end
%% Display
cmap = jet(25);

figure('Color', 'w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile

% subplot(221)
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
% subplot(222)
pcolor(D.lon, D.lat, varRaw_s')
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
% subplot(223)
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
% subplot(224)
pcolor(D.lon, D.lat, varRaw_b')
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




