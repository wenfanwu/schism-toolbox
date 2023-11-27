function check_schism_init(Mobj, DS, InitCnd, varName, iDep)
% Check the interpolation for the initial inputs

%%  Parse inputs
if nargin < 4
    varName = 'temp';
end
if nargin < 5
    iDep = 1;
end

D = DS.(varName);
switch varName
    case 'ssh'
        varRaw_b = D.var;
        varNew_b = InitCnd.ssh;
        varRaw_s = D.var;
        varNew_s = InitCnd.ssh;
    otherwise
        msk3d = isnan(Mobj.depLayers)';
        InitCnd.(varName)(msk3d) = nan;

        varRaw = D.var;
        varNew = InitCnd.(varName);

        varRaw2d = reshape(varRaw, [size(varRaw, 1)*size(varRaw,2) size(varRaw,3)]);
        ind_btm2d = sum(~isnan(varRaw2d'));
        ind_btm2d(ind_btm2d==0) = 1;
        varBtm_raw2d = arrayfun(@(x) varRaw2d(x, ind_btm2d(x)), 1:size(varRaw2d,1));

        varRaw_b = reshape(varBtm_raw2d, size(varRaw, [1 2]));
        varRaw_s = squeeze(varRaw(:,:, min(iDep, numel(D.depth))));

        ind_btm = sum(~isnan(Mobj.depLayers));
        varNew_b = arrayfun(@(x) varNew(x, ind_btm(x)), 1:Mobj.nNodes);
        varNew_s = varNew(:, min(iDep, Mobj.maxLev));
end
%% Display
cmap = jet(25);

figure('Color', 'w')
% tiledlayout(2,2,'TileSpacing','tight')  % better alternatives for advanced MATLAB version
% nexttile

subplot(221)
disp_schism_var(Mobj, varNew_s)
hold on
plot_schism_bnds(Mobj, [1 1], 'Color', 'k')
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
colorbar
colormap(cmap)
varLim = caxis;
title('SCHISM (surface)', 'FontWeight','bold')

% nexttile
subplot(222)
pcolor(D.lon, D.lat, varRaw_s')
shading flat
hold on
plot_schism_bnds(Mobj, [1 1], 'Color', 'k')
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
colorbar
colormap(cmap)
caxis(varLim)
xlabel('Longitude (°E)', 'FontWeight','bold')
ylabel('Latitude (°N)', 'FontWeight','bold')
title('Raw Data (surface)', 'FontWeight','bold')

% nexttile
subplot(223)
disp_schism_var(Mobj, varNew_b)
hold on
plot_schism_bnds(Mobj, [1 1], 'Color', 'k')
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
colorbar
colormap(cmap)
varLim = caxis;
title('SCHISM (bottom)', 'FontWeight','bold')

% nexttile
subplot(224)
pcolor(D.lon, D.lat, varRaw_b')
shading flat
hold on
plot_schism_bnds(Mobj, [1 1], 'Color', 'k')
axis image
box on
xlim(Mobj.region(1:2))
ylim(Mobj.region(3:4))
colorbar
colormap(cmap)
caxis(varLim)
xlabel('Longitude (°E)', 'FontWeight','bold')
ylabel('Latitude (°N)', 'FontWeight','bold')
title('Raw Data (bottom)', 'FontWeight','bold')
% auto_center

end




