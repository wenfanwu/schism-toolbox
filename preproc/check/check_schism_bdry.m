function check_schism_bdry(Mobj, DS, BdryCnd, varName, iTime)
% Check the interpolation at the boundary nodes
% 
%% Syntax
% 
%
%% Description 
% 
%
%% Examples
%
%
%% Input Arguments
%
%
%% Output Arguments
% 
% 
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-10-22.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if nargin < 4
    varName = 'temp';
end
if nargin < 5
    iTime = 1;
end
D = DS.(varName);

nNodes_obc = Mobj.nNodes_obc;
nDeps_new = Mobj.maxLev;
nDeps_raw = numel(D.depth);

varData = BdryCnd.(varName);
varData2 = D.var;

varNew = squeeze(varData(:,:,iTime))';
depNew = Mobj.depLayers(:, Mobj.obc_nodes_tot);
distNew = repmat(1:nNodes_obc, nDeps_new, 1);

indTime = minfind(D.time, Mobj.time(iTime));
varRaw = squeeze(varData2(:,:, indTime))';
depRaw = repmat(D.depth(:), 1, nNodes_obc);
distRaw =  repmat(1:nNodes_obc, nDeps_raw, 1);
%% Display
depRaw = -abs(depRaw);
depNew = -abs(depNew);

figure('Color', 'w')
subplot(211)
pcolor(distNew, depNew, varNew)
shading flat
colormap(jet(25))
colorbar
varLim = caxis;
yvarLim = ylim;
xlabel('Along Open Boundary Nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
title(['SCHISM (', datestr(Mobj.time(iTime), 'yyyy-mm-dd'), ')'])

subplot(212)
pcolor(distRaw, depRaw, varRaw)
shading flat
colormap(jet(25))
colorbar
caxis(varLim)
ylim(yvarLim)
xlabel('Along Open Boundary Nodes', 'FontWeight','bold')
ylabel('Depth (m)', 'FontWeight','bold')
title(['Raw Data (', datestr(D.time(indTime), 'yyyy-mm-dd'), ')'])

end













