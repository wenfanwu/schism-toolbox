function angles = calc_schism_angles(Mobj)
% Calculate the interior angles (deg) of each cell
%
%% Syntax
% angles = calc_schism_angles(Mobj)
%
%% Description
% angles = calc_schism_angles(Mobj) calculates the interior angles of each cell
%
%% Examples
% angles = calc_schism_angles(Mobj);
% angles_sum = sum(angles, 2, 'omitnan');
% 
% figure('Position', [680, 182, 592, 696], 'Color', 'w')
% disp_schism_var(Mobj, angles_sum, 'EdgeColor', 'k')
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
%
%% Output Arguments
% angles - interior angles; double
%       the interior angles (0-180) of each cell.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 24 Nov 2024.
% Email: wwu@vims.edu
%
% See also: calc_schism_skew

%% Parse inputs
angles = nan(Mobj.nElems,4);

tri = Mobj.tri;
tri(Mobj.i34==3,4) = tri(Mobj.i34==3,3);  % fillmissing

lon_tri = Mobj.lon(tri);
lat_tri = Mobj.lat(tri);

% vertex #1
v1 = [lon_tri(:,4), lat_tri(:,4)] - [lon_tri(:,1), lat_tri(:,1)];
v2 = [lon_tri(:,2), lat_tri(:,2)] - [lon_tri(:,1), lat_tri(:,1)];
angles(:, 1) = calc_vector_angles(v1, v2);

% vertex #2
v1 = [lon_tri(:,1), lat_tri(:,1)] - [lon_tri(:,2), lat_tri(:,2)];
v2 = [lon_tri(:,3), lat_tri(:,3)] - [lon_tri(:,2), lat_tri(:,2)];
angles(:, 2) = calc_vector_angles(v1, v2);

% vertex #3
v1 = [lon_tri(:,2), lat_tri(:,2)] - [lon_tri(:,3), lat_tri(:,3)];
v2 = [lon_tri(:,4), lat_tri(:,4)] - [lon_tri(:,3), lat_tri(:,3)];
angles(:, 3) = calc_vector_angles(v1, v2);

% vertex #4
v1 = [lon_tri(:,3), lat_tri(:,3)] - [lon_tri(:,4), lat_tri(:,4)];
v2 = [lon_tri(:,1), lat_tri(:,1)] - [lon_tri(:,4), lat_tri(:,4)];
angles(:, 4) = calc_vector_angles(v1, v2);

% adjust for triangular cells
angles(Mobj.i34==3,3) = 180-angles(Mobj.i34==3,1)-angles(Mobj.i34==3,2);
angles(Mobj.i34==3,4) = nan;

end

function A = calc_vector_angles(v1, v2)
% Calculate the angles (deg) between two vectors

% dot products and vector modules
dot_prot = dot(v1, v2, 2);
norm_v1 = sqrt(sum(v1.^2, 2));
norm_v2 = sqrt(sum(v2.^2, 2));

% the angle between vectors (deg)
cos_theta = dot_prot./(norm_v1.*norm_v2);
cos_theta = round(cos_theta,6);
cos_theta(isnan(cos_theta)) = 1;

A = acos(max(-1, min(1, cos_theta)))*180/pi;
end

