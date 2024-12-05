function skn = calc_schism_skew(Mobj, mtype)
% Calculate the skewness (0-1) of each cell
%
%% Syntax
% skn = calc_schism_skew(Mobj)
%
%% Description
% skn = calc_schism_skew(Mobj) calculates the skewness of each cell
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
% mtype - method type; double
%       the method type determing the skewness definition method. Default:
%       mtype = 1; Three methods are provided for skewness calculation. 
%       1) the deviation to equilateral-triangle/rectangle
%               skn ranges in 0-1, with lower values meaning higher quality
%       2) the ratio between the min and max interior angles
%               skn ranges in 0-1, with higher values meaning higher quality
%       3) the ratio between the max. side and the equivalent radius
%               skn is largerr than 1, with lower values meaning higher quality
%
%% Output Arguments
% skn - skewness; double
%       the skewness of each cell. 
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 24 Nov 2024.
% Email: wwu@vims.edu
%
% See also: calc_schism_angles

%% Parse inputs
if nargin < 2
    mtype = 1;
end

switch mtype
    case 1  % the deviation from equilateral-triangle or rectangle.
        angles = calc_schism_angles(Mobj);
        i34 = Mobj.i34;
        skn = nan(Mobj.nElems,1);
        skn(i34==3) = sqrt(sum((angles(i34==3)/60-1).^2, 2, 'omitnan')./3);
        skn(i34==4) = sqrt(sum((angles(i34==4)/90-1).^2, 2, 'omitnan')./4);
        % normalized to 0-1 (by dividing the theoretical max. value)
        skn(i34==3) = skn(i34==3)/sqrt(2);
        skn(i34==4) = skn(i34==4)/sqrt(3);

    case 2  % the ratio between the min and max interior angles
        angles = calc_schism_angles(Mobj);
        max_angle = max(angles, [], 2, 'omitnan');
        min_angle = min(angles, [], 2, 'omitnan');
        skn = min_angle./max_angle;

    case 3  % the ratio between the max. side and the equivalent radius
        edges = calc_schism_edge(Mobj);
        R = calc_schism_reso(Mobj, 4); % circumcircle radius
        max_edge = max(edges, [], 2, 'omitnan');
        skn = max_edge./R;
        
end

end