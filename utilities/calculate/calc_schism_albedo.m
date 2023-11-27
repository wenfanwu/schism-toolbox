function albedo_val = calc_schism_albedo(Mobj, ntype, disp_flag)
% Calculate the albedo based on different empirical formulae
% 
%% Syntax
% albedo_val = calc_schism_albedo(Mobj, ntype)
% albedo_val = calc_schism_albedo(Mobj, ntype, disp_flag)
% 
%% Description 
% albedo_val = calc_schism_albedo(Mobj, ntype) calculates the albedo
% based on empirical formula.
% albedo_val = calc_schism_albedo(Mobj, ntype, disp_flag) displays the albedo
% on the unstructured mesh.
%
%% Examples
% albedo_val = calc_schism_albedo(Mobj, 1, 'on');
%
%% Input Arguments
% Mobj --- the mesh object
% ntype --- index of selected formula
% disp_flag --- display the calculated results or not ('on' or 'off')
%
%% Output Arguments
% albedo_val --- the calculated albedo values
% 
%% Notes
% Users can easily add your own empirical formula in the code 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-23.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if nargin<2
    ntype = 1;
end
if nargin<3
    disp_flag = 'off';
end

switch ntype
    case 1  % Large and Yeager (2009)
        albedo_val = 0.069-0.011*cosd(2*Mobj.lat); 
    otherwise
        % add your own formula here
end

%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, albedo_val)
    colormap(jet(25))
    box on
    axis image
end
end

% Reference:
% Large, W., & Yeager, S. G. (2009). The global climatology of an
% interannually varying air–sea flux data set. Climate dynamics, 33(2),
% 341-364.