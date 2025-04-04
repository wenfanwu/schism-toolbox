function check_schism_bdry(Mobj, DS, BdryCnd, varName, idx_time, obc_bnds)
% Check the interpolation at the boundary nodes
%
%% Syntax
% check_schism_bdry(Mobj, DS, BdryCnd, varName)
% check_schism_bdry(Mobj, DS, BdryCnd, varName, idx_time)
% check_schism_bdry(Mobj, DS, BdryCnd, varName, idx_time, obc_bnds)
%
%% Description
% check_schism_bdry(Mobj, DS, BdryCnd, varName) checks the interpolation along open boundary nodes
% check_schism_bdry(Mobj, DS, BdryCnd, varName, idx_time) specifies the time steps
% check_schism_bdry(Mobj, DS, BdryCnd, varName, idx_time, obc_bnds) specifies the open boundaries
%
%% Examples 
% DS = prep_schism_bdry(Mobj, 'hycom_bys'); 
% varList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'}; 
% BdryCnd = interp_schism_bdry(Mobj, DS, varList);
% 
% check_schism_bdry(Mobj, DS, BdryCnd, 'temp',1)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% DS - raw data struct; datastruct
%       the datastruct containing boundary data, resulting from "prep_schism_bdry"
% BdryCnd - data object; datastruct
%       the datastruct containing boundary data, resulting from "interp_schism_bdry"
% idx_time - index time (optional); numeric
%       the index of time steps. Default: idx_time = 1.
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
% Last Updated on 1 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: check_schism_init

%% Parse inputs
if nargin < 4; varName = 'temp'; end
if nargin < 5; idx_time = 1; end

D = DS.(varName);
% Determine the used open boundary segments
obc_bnds_all = find(ismember(Mobj.obc_nodes(1,:), D.ind));

if nargin < 6; obc_bnds = 'all'; end
if strcmpi(obc_bnds, 'all'); obc_bnds = obc_bnds_all; end

idx_bnds = ismember(obc_bnds, obc_bnds_all);
if numel(find(idx_bnds==0))>0
    error('some open boundary segments in obc_bnds are unused or inexistent!')
end
%% Prepare data
obc_nodes = Mobj.obc_nodes(:, obc_bnds); 
obc_nodes = obc_nodes(:);
obc_nodes(obc_nodes==0) = [];
nps = numel(obc_nodes);

nz_new = Mobj.maxLev;
nz_raw = numel(D.depth);

v1 = BdryCnd.(varName);
v2 = D.var; 

var_new = squeeze(v1(:,:,idx_time))';
dep_new = -abs(Mobj.depLayers(:, obc_nodes));
dist_new = repmat(1:nps, nz_new, 1);

idx_time2 = minfind(D.time, Mobj.time(idx_time));  % time index in the raw data.
var_raw = squeeze(v2(:,:, idx_time2))';
dep_raw = -abs(repmat(D.depth(:), 1, nps));
dist_raw =  repmat(1:nps, nz_raw, 1);

%% Display
figure('Color', 'w')
subplot(211)
pcolor(dist_new, dep_new, var_new)
shading flat
colormap(jet(25))
colorbar
cm = caxis; %#ok<*CAXIS>
ym = ylim;
xlabel('Along open boundary nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
set(gca, 'Layer', 'top')
title(['SCHISM (', datestr(Mobj.time(idx_time), 'yyyy-mm-dd'), ')']) %#ok<*DATST>

subplot(212)
pcolor(dist_raw, dep_raw, var_raw)
shading flat
colormap(jet(25))
colorbar
caxis(cm)
ylim(ym)
xlabel('Along open boundary nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
set(gca, 'Layer', 'top')
title(['Raw Data (', datestr(D.time(idx_time2), 'yyyy-mm-dd'), ')'])

end













