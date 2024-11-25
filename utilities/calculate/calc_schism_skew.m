function skn = calc_schism_skew(Mobj)
% Calculate the skewness (0-1) of each cell
%
%% Syntax
% skn = calc_schism_angles(Mobj)
%
%% Description
% skn = calc_schism_angles(Mobj) calculates the skewness of each cell
%
%% Examples
% skn = calc_schism_skew(Mobj);
% skn(skn<0.25) = nan;
% 
% figure('Position', [447, 428, 736, 464], 'Color', 'w')
% disp_schism_var(Mobj, skn, 'EdgeColor', 'k')
% axis image
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
%
%% Output Arguments
% skn - skewness; double
%       the skewness (0-1) of each cell. For triangular/quadrangular cell,
%       skewness is defined as its deviation toequilateral-triangle/rectangle, 
%       skn = 0 for equilateral-triangle or rectangle.  
%       the formula is given below: 
%       skn_tri = sqrt(((A/60-1)^2 + (B/60-1)^2 + (C/60-1)^2)/3)/sqrt(2)
%       skn_quad = sqrt(((A/90-1)^2 + (B/90-1)^2 + (C/90-1)^2 + (D/90-1)^2)/4)/(sqrt(3))
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 24 Nov 2024.
% Email: wwu@vims.edu
%
% See also: calc_schism_angles

%% Parse inputs
angles = calc_schism_angles(Mobj);

i34 = Mobj.i34;
skn = nan(Mobj.nElems,1);
skn(i34==3) = sqrt(sum((angles(i34==3)/60-1).^2, 2, 'omitnan')./3);
skn(i34==4) = sqrt(sum((angles(i34==4)/90-1).^2, 2, 'omitnan')./4);

% normalized to 0-1 (by dividing the theoretical max. value)
skn(i34==3) = skn(i34==3)/sqrt(2); 
skn(i34==4) = skn(i34==4)/sqrt(3); 

end