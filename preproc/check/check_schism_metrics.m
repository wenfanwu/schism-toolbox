function check_schism_metrics(Mobj, CFL_limit)
% Check the horizontal grid quality
% 
%% Syntax
% check_schism_metrics(Mobj)
% check_schism_metrics(Mobj, 0.4)
% 
%% Description 
% check_schism_metrics(Mobj) check the CFL restriction.
% 
% check_schism_metrics(Mobj) set the minimum CFL number on
% your own.
% 
%% Input Arguments
% Mobj --- mesh object
% CFL_val --- min. CFL numbers
%
%% Output Arguments
% None
% 
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2023-11-26.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: calc_schism_CFL

%% Parse inputs
if nargin < 2
    CFL_limit = 0.4;
end
g = 9.80665;
ua = 0;
h = abs(Mobj.depth);
dt = Mobj.dt;

%% Calculation
R = calc_schism_cradius(Mobj);      % use the circumradius
% R = calc_schism_sidelen(Mobj);  % use the side length

dx = R(:);
CFL_val = (ua+sqrt(g*h)).*dt./dx;

dx0 = sqrt(g*h).*dt./CFL_limit;

dx_diff = dx-dx0;
dx_diff(dx_diff<=0) = nan;

%% Display
cmap = jet(15);

figure('Color', 'w')
% tiledlayout(2,2,'TileSpacing','tight')   % better alternative for advanced versions of Matlab

% nexttile
subplot(221)
disp_schism_var(Mobj, dx/1000)
axis image
box on
title('Actual grid resolutions (km)')

% nexttile
subplot(222)
disp_schism_var(Mobj, dx0/1000)
axis image
box on
title('Acceptable coarsest grid resolutions(km)')

% nexttile
subplot(223)
disp_schism_var(Mobj, CFL_val)
axis image
box on
caxis([0 1])
title('CFL numbers')

% nexttile
subplot(224)
disp_schism_var(Mobj, dx_diff/1000, 'EdgeColor', 'k')
axis image
box on
colormap(cmap)
title('Excessive resolutions (km)')

auto_center
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


