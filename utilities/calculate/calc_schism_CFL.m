function calc_schism_CFL(Mobj, CFL_limit)
% Check the inverse CFL constraints.
%
%% Syntax
% calc_schism_CFL(Mobj, CFL_limit)
% 
%% Description 
% calc_schism_CFL(Mobj, CFL_limit) returns a graph showing the relationship
%       between horizontal resolution and water depth.
% calc_schism_CFL(Mobj, CFL_limit) sets the critical CFL number
% 
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% CFL_limit - CFL limit; double
%       the critical CFL number. by default: CFL_limit = 0.4;
% 
%% Notes
% the returned diagram is just a reference for the grid quality.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 11 Nov 2024.
% Email: wwu@vims.edu
% 
% See also: check_schism_metrics

%% Parse inputs
if nargin==1
    CFL_limit = 0.4;
end

uflow = 0;
g = 9.80665;
dt = Mobj.dt;
depthc = max(Mobj.depthc,0);  % avoid negative depth (e.g., hills).

calc_dx = @(h,dt) sqrt(uflow+g*h).*dt./CFL_limit;

h0 = min(depthc(:)):0.25:max(depthc(:));
dx0 = calc_dx(h0, dt);

R = calc_schism_reso(Mobj);
Rc = calc_dx(depthc, dt);
%% Display
figure('Color', 'w')
plot(h0, dx0/1000, 'LineWidth', 1, 'Color', 'k')
box on; grid on
xlabel('Depth (m)', 'FontWeight', 'bold')
ylabel('Horizontal resolution (km)', 'FontWeight', 'bold')
hold on
scatter(depthc(R>Rc), R(R>Rc)/1e3, 2, 'red', 'filled')
scatter(depthc(R<=Rc), R(R<=Rc)/1e3, 2, 'blue', 'filled')

rb = numel(find(R>Rc))/Mobj.nElems*100;
rg = numel(find(R<=Rc))/Mobj.nElems*100;

title(['dt = ', num2str(dt, '%02d'), 's'])
legend({'theoretical coarsest resolution', ['unqualified cells (', num2str(rb,'%.02f'),'%)'],  ['qualified cells (', num2str(rg,'%.02f'),'%)']})
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

