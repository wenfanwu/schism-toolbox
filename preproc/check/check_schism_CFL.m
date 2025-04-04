function check_schism_CFL(Mobj, CFL_limit)
% Check the inverse CFL constraints for SCHISM grid.
% 
%% Syntax
% check_schism_CFL(Mobj)
% check_schism_CFL(Mobj, CFL_limit)
% 
%% Description 
% check_schism_CFL(Mobj)checks the inverse CFL constraints
% check_schism_CFL(Mobj, CFL_limit) specifies the critical CFL number.
% 
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct containing mesh info.
% CFL_limit - the critical CFL number; numeric
%       the critical CFL number. Default: CFL_limit = 0.4
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 3 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: check_schism_hydrostatic

%% Parse inputs
if nargin < 2; CFL_limit = 0.4; end

% define essential parameters
g = 9.80665;
uvel = 0;
h = max(Mobj.depthc(:), 0.1);   % ignore extremely shallow areas
dt = Mobj.dt;
dx = calc_schism_reso(Mobj); % horizontal resolution (m)

%% Calculate CFL numbers
CFL_val = (uvel+sqrt(g*h)).*dt./dx(:);

% the theoratical coarest dx at given dt and CFL number.
calc_dx = @(h, dt) (uvel+sqrt(g*h))*dt./CFL_limit; 

h_s = min(h(:)):0.25:max(h(:));

dx_s = calc_dx(h_s, dt);
dx_c = calc_dx(h, dt);

% the over-refined zones
dx_diff = dx(:) - dx_c(:);
dx_diff(dx_diff<=0) = nan;

%% Figure-1: the spatial map of inverse CFL constraints.
figure('Color', 'w')
tiledlayout(2,2,'TileSpacing','tight')   % better alternative of subplot since R2019b

nexttile
% subplot(221)
disp_schism_var(Mobj, dx/1000)
axis image
box on
title('Horizontal resolutions (km)')

nexttile
% subplot(222)
disp_schism_var(Mobj, dx_c/1000)
axis image
box on
title('Theoreticel coarsest resolutions (km)')

nexttile
% subplot(223)
disp_schism_var(Mobj, CFL_val)
axis image
box on
caxis([0 1]) %#ok<*CAXIS>
title(['CFL numbers', ' (CFL > ', num2str(CFL_limit), ')'])

nexttile
% subplot(224)
disp_schism_var(Mobj, dx_diff/1000, 'EdgeColor', 'k')
axis image
box on
colormap(jet(15))
title('Excessive resolutions (km)')

%% Figure-2: the theoretical coarsest resolutions as a function of water depth.
figure('Color', 'w')
plot(h_s, dx_s/1000, 'LineWidth', 1, 'Color', 'k')
box on; grid on
xlabel('Depth (m)', 'FontWeight', 'bold')
ylabel('Horizontal resolution (km)', 'FontWeight', 'bold')
hold on
scatter(h(dx>dx_c), dx(dx>dx_c)/1e3, 2, 'red', 'filled')
scatter(h(dx<=dx_c), dx(dx<=dx_c)/1e3, 2, 'blue', 'filled')
auto_center

r1 = numel(find(dx>dx_c))/Mobj.nElems*100;
r2 = numel(find(dx<=dx_c))/Mobj.nElems*100;

title(['dt = ', num2str(dt, '%02d'), 's'])
legend({'theoretical coarsest resolution', ['unqualified cells (', num2str(r1,'%.02f'),'%)'],  ['qualified cells (', num2str(r2,'%.02f'),'%)']})

end

% Table 5.1: Coarsest grid resolution at sample depths, assuming a ‘worse
% case’ scenario of t=100s.
%   h (m)     xmax (m)
%   1             790
%   10           2500
%   50           5.5e3
%   100         7.9e3
%   500         1.7e4
%   1000        2.5e4
%   4000        5e4


