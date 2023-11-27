function Cd = calc_schism_bfric(Mobj, ntype, tuning_coefs, disp_flag)
% Calculate the bottom friction using different empirical formulae
%
%% Syntax
% Cd = calc_schism_bfric(Mobj, ntype, tuning_coefs)
%
%% Description
% Cd = calc_schism_bfric(Mobj, ntype, tuning_coefs) returns the Cd (bottom
% drag coefficient).
%
%% Input Arguments
% Mobj --- the mesh object
% ntype --- the formula type to calculate the Cd.
% turning_coefs --- turning coefficient pair. tuning_coefs = [fmc, h_min];
% fmc is the factor of maning coefficient; h_min is the minmum depth in
% manning formula.
%
%% Output Arguments
% Cd --- bottom drag coefficient.
%
%% Notes
% If neccessary, add your own empirical formula in the code
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2 Dec. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: disp_schism_var

%% Parse inputs
if nargin < 2
    ntype = 1;
end
if nargin < 3
    tuning_coefs = [0.025 3];
end
if nargin < 4
    disp_flag = 'on';
end

%% Caculate
switch ntype
    case 1  % Method-1: manning method
        fmc = tuning_coefs(1);
        h_min = tuning_coefs(2);
        fcn = @(x) (fmc^2*9.810)./(max(x, h_min).^(1.0/3.0));
    case 2
        fmc = tuning_coefs(1);
        h_min = tuning_coefs(2);
        fcn = @(x) (fmc^2*9.810)./(max(x, h_min).^(2.0/3.0));
end
Cd = fcn(Mobj.depth);

%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w');
    disp_schism_var(Mobj, Cd)
    colormap(jet(25))
    axis image
    box on
end
end


