function check_schism_bdry(Mobj, DS, BdryCnd, varName, timeTick, obc_bnds)
% Check the interpolation at the boundary nodes
%
%% Syntax
% check_schism_bdry(Mobj, DS, BdryCnd, varName)
% check_schism_bdry(Mobj, DS, BdryCnd, varName, timeTick)
% check_schism_bdry(Mobj, DS, BdryCnd, varName, timeTick,  obc_bnds)
%
%% Description
% check_schism_bdry(Mobj, DS, BdryCnd, varName) checks the interpolation along open boundary nodes
% check_schism_bdry(Mobj, DS, BdryCnd, varName, timeTick) specifies the time moment
% check_schism_bdry(Mobj, DS, BdryCnd, varName, timeTick, obc_bnds) specifies the open boundaries
%
%% Examples 
% DS = prep_schism_bdry(Mobj, 'hycom_bys'); 
% varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; 
% BdryCnd = interp_schism_bdry(Mobj, DS, varList);
% 
% check_schism_bdry(Mobj, DS, BdryCnd, 'temp', Mobj.time(1))
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% DS - raw data struct; datastruct
%       the datastruct containing boundary data, resulting from "prep_schism_bdry"
% BdryCnd - data object; datastruct
%       the datastruct containing boundary data, resulting from "interp_schism_bdry"
% timeTick - time moment (optional); datetime
%       the time moment to show the results. Default: timeTick = Mobj.time(1)
% obc_bnds - open boundary index (optional); numeric
%       the index of open boundary segments. e.g., obc_bnds =
%       1:Mobj.obc_counts. All available open boundary segments in "DS"
%       will be used by default.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022. 
% Last Updated on 18 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: check_schism_init

%% Parse inputs
if nargin < 4; varName = 'temp'; end
if nargin < 5; timeTick = Mobj.time(1); end

D1 = DS(strcmp({DS.Variable}, varName));
D2 = BdryCnd(strcmp({BdryCnd.Variable}, varName));

%% Determine the open boundaries
obc_bnds_all = find(ismember(Mobj.obc_nodes(1,:), D1.Nodes));
if nargin < 6; obc_bnds = obc_bnds_all; end
idx_bnds = ismember(obc_bnds, obc_bnds_all);
if numel(find(idx_bnds==0))>0
    error('some open boundary segments in obc_bnds are unused or inexistent!')
end

%% Prepare data
obc_nodes = Mobj.obc_nodes(:, obc_bnds);  obc_nodes = obc_nodes(:); obc_nodes(obc_nodes==0) = [];
nps = numel(obc_nodes);

nz_new = Mobj.maxLev; nz_raw = numel(D1.Depth);
v1 = D1.Data; v2 = D2.Data; 

idx_time1 = minfind(D2.Time, timeTick);
var_new = squeeze(v2(:,:,idx_time1));
dep_new = fillmissing(-abs(Mobj.depLayers(:, obc_nodes)), 'previous',1);
dist_new = repmat(1:nps, nz_new, 1);

idx_time2 = minfind(D1.Time, Mobj.time(idx_time1));  % time index in the raw data.
var_raw = squeeze(v1(:,:, idx_time2));
dep_raw = -abs(repmat(D1.Depth(:), 1, nps));
dist_raw =  repmat(1:nps, nz_raw, 1);

%% Display
bdep = -abs(Mobj.depth(obc_nodes));

figure('Color', 'w')
subplot(211)
pcolor(dist_new, dep_new, var_new)
shading flat
hold on
plot(bdep, 'LineWidth', 1, 'Color', 'k')
colormap(jet(25))
colorbar
cm = caxis; %#ok<*CAXIS>
ym = ylim;
xlabel('Along open boundary nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
set(gca, 'Layer', 'top')
title(['SCHISM (', datestr(Mobj.time(idx_time1), 'yyyy-mm-dd HH:MM:SS'), ')']) %#ok<*DATST>

is_nan = dep_raw<repmat(bdep(:)', [size(dep_raw,1) 1]);
var_raw(is_nan) = nan;

subplot(212)
pcolor(dist_raw, dep_raw, var_raw)
shading flat
hold on
plot(bdep, 'LineWidth', 1, 'Color', 'k')
colormap(jet(25))
colorbar
caxis(cm)
ylim(ym)
xlabel('Along open boundary nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
set(gca, 'Layer', 'top')
title(['Raw Data (', datestr(D1.Time(idx_time2), 'yyyy-mm-dd HH:MM:SS'), ')'])

end













