function albedo_val = calc_schism_albedo(Mobj, mtype, disp_flag)
% Calculate the albedo based on different empirical formulae
% 
%% Syntax
% albedo_val = calc_schism_albedo(Mobj, mtype)
% albedo_val = calc_schism_albedo(Mobj, mtype, disp_flag)
% 
%% Description 
% albedo_val = calc_schism_albedo(Mobj, mtype) calculates the albedo based
%       on empirical formula. 
% albedo_val = calc_schism_albedo(Mobj, mtype, disp_flag) displays the
%       albedo or not.
%
%% Examples
% albedo_val = calc_schism_albedo(Mobj, 1, 'on');
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct used to store the mesh info.
% mtype - the method type; numeric
%       the method used to calculate albedo; Default: mtype = 1;
%       mtype = 1 (Large and Yeager, 2009): calculate the albedo as a
%       function of latitude;
% disp_flag - display flag; char
%       the flat used to display the resulting albedo or not ('on' or
%       'off'); Default: disp_flag = 'off';
%
%% Output Arguments
% albedo_val - albedo values; numeric
%       the calculated albedo values.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 28 May 2025.
% Email: wwu@vims.edu
% 
% See also: calc_schism_bfric

%% Parse inputs
if nargin<2; mtype = 1; end
if nargin<3; disp_flag = 'off'; end

switch mtype
    case 1  % Large and Yeager (2009)
        albedo_val = 0.069-0.011*cosd(2*Mobj.lat); 
    otherwise
        % add your own formula here
end

%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, albedo_val)
    hold on
    plot_schism_bnds(Mobj)
    colormap(jet(25))
    box on
    axis image
    auto_center
end
end

% Reference:
% Large, W., & Yeager, S. G. (2009). The global climatology of an
% interannually varying air–sea flux data set. Climate dynamics, 33(2),
% 341-364.