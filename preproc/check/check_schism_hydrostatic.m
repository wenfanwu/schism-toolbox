function check_schism_hydrostatic(Mobj, aspect_ratio, disp_flag)
% Check the hydrostatic assumption 
% 
%% Syntax
% check_schism_hydrostatic(Mobj)
% check_schism_hydrostatic(Mobj, aspect_ratio)
% check_schism_hydrostatic(Mobj, aspect_ratio, disp_flag)
% 
%% Description 
% check_schism_hydrostatic(Mobj) checks the hydrostatic assumption
% check_schism_hydrostatic(Mobj, aspect_ratio) specifies the aspect ratio (L/H)
% check_schism_hydrostatic(Mobj, aspect_ratio, disp_flag) displays the results or not.
%
%% Examples
% check_schism_hydrostatic(Mobj, 10, 'on')
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct containing mesh info.
% aspect_ratio - the critical aspect ratio; numeric
%       the critical aspect ratio between horizontal resolution and local
%       water depth.
% disp_flag - display flag; char
%       the flag that determines to display the results or not ('on'/'off'); 
%       Default: disp_flag = 'on'.
%
%% Output Arguments
% None
% 
%% Notes
% Make sure the horizontal scale >> vertical scale, otherwise the
% hydrostatic assumption will be violated and there might be spurious
% upwelling.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 3 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: check_schism_CFL

%% Parse inputs
if nargin < 2
   aspect_ratio = 10;  % Note that this aspect ratio is quite strict actually.
end
if nargin < 3
   disp_flag = 'on';
end

%% Calculation
R = calc_schism_reso(Mobj);   

h = max(Mobj.depthc(:), 0.1);
Rd = R(:)./h; Rd(Rd>aspect_ratio) = nan;

N = numel(find(~isnan(Rd)));
%% Display
switch lower(disp_flag)
    case 'on'
        figure('Color', 'w')
        disp_schism_var(Mobj, Rd, 'EdgeColor', 'k')
        axis image
        box on
        caxis([0 aspect_ratio]) %#ok<CAXIS>
        title(['# of violating nodes = ', num2str(N), ' (L/H < ', num2str(aspect_ratio), ')'])
        auto_center
end

end