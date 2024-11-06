function Rd = check_schism_hydrostatic(Mobj, cut_val, disp_flag)
% Check the hydrostatic assumption 
% 
%% Syntax
% Rd = check_schism_hydrostatic(Mobj)
% Rd = check_schism_hydrostatic(Mobj, cut_val)
% Rd = check_schism_hydrostatic(Mobj, cut_val, disp_flag)
% 
%% Description 
% Rd = check_schism_hydrostatic(Mobj) checks the hydrostatic assumption
% Rd = check_schism_hydrostatic(Mobj, cut_val) specifies the critical ratio
% Rd = check_schism_hydrostatic(Mobj, cut_val, disp_flag) displays the
% results or not.
%
%% Examples
% Rd = check_schism_hydrostatic(Mobj, 10, 'on')
%
%% Input Arguments
% Mobj --- the mesh object
% cut_val --- the critical ratio between horizontal resolution and local
% water depth
% disp_flag --- display the results or not ('on'/'off')
%
%% Output Arguments
% Rd --- index of notes that violate the hydrostatic assumption
% 
%% Notes
% Make sure the horizontal scale >> vertical scale, otherwise the
% hydrostatic assumption will be violated and there might be spurious
% upwelling. By default, it holds when horizontal resolution is 10 times
% larger than local water depth.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-23.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 
%% Parse inputs
if nargin < 2
   cut_val = 10;
end
if nargin < 3
   disp_flag = 'on';
end

%% Calculation
R = calc_schism_reso(Mobj);   

Rd = R(:)./abs(Mobj.depthc(:));
Rd(Rd>cut_val) = nan;

N = numel(find(~isnan(Rd)));

%% Display
switch disp_flag
    case 'on'
        figure('Color', 'w')
        disp_schism_var(Mobj, Rd, 'EdgeColor', 'k')
        hold on
        plot_schism_bnds(Mobj)
        axis image
        box on
        caxis([0 cut_val]) %#ok<CAXIS>
        title(['# of violating points = ', num2str(N)])
        auto_center
end

end