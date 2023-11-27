function calc_schism_CFL(Mobj, CFL_limit)
% Estimate the inversed CFL constraints
%
%% Syntax
% calc_schism_CFL(Mobj, CFL_val)
% 
%% Description 
% calc_schism_CFL(Mobj, CFL_val) returns a graph showing the relationship
% between horizontal resolution and water depth.
% calc_schism_CFL(Mobj, CFL_val) sets the critical CFL number
% 
%% Input Arguments
% Mobj --- mesh object
% CFL_val --- critical CFL number
%
%% Output Arguments
% None
% 
%% Notes
% this function needs improvements since it is time-consuming!
% 
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: check_schism_metrics
%
%% Parse inputs
if nargin==1
    CFL_limit = 0.4;
end

uflow = 0;
g = 9.8;
dt = Mobj.dt;

fcn = @(h,dt) sqrt(uflow+g*h).*dt./CFL_limit;

h0 = min(Mobj.depth(:)):max(Mobj.depth(:));
dx0 = fcn(h0, dt);

%% Display
figure('Color', 'w')
plot(h0, dx0/1000, 'LineWidth', 1, 'Color', 'b')
box on; grid on
xlabel('Depth (m)', 'FontWeight', 'bold')
ylabel('Max. grid size (km)', 'FontWeight', 'bold')

end



