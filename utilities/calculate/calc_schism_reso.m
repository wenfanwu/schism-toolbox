function R = calc_schism_reso(Mobj, mtype)
% Calculate horizontal resolution of grid cells
%
%% Syntax
% R = calc_schism_reso(Mobj)
% R = calc_schism_reso(Mobj, mtype)
%
%% Description
% R = calc_schism_reso(Mobj) calcuates the grid resolutions.
% R = calc_schism_reso(Mobj, mtype) specifies the method to calculate grid resolutions.
%
%% Example
% R = calc_schism_reso(Mobj, 1)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% mtype - method type; numeric
%       mtype determines the definition method of grid resolution.
%       Five methods are provided (default: mtype=4):
%       1) average side length           (x1+x2+x3)/3 or (x1+x2+x3+x4)/4
%       2) equal area circle radius     sqrt(S/pi)
%       3) square root area                sqrt(S)
%       4) circumcircle radius             r1=(x1*x2*x3)/(4S);
%       5) incircle radius                     r2=2S/(x1+x2+x3)
%
%% Output Arguments
% R - horizontal resolution; numeric
%       the horizontal/spatial resolution of each grid cell.
%
%% Notes
% All quadrangular cells will be splited into two triangular cells to calculate
% the radius independently. As for quadrangular cells, the circumcircle
% radius is the largest one from its two triangular cells. and the incircle
% radius is the sum of its two triangular cells' incircle radius. 
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 11 Nov 2024.
% Email: wwu@vims.edu
%
% See also: calc_schism_area and calc_schism_edge

%% Parse inputs
if nargin == 1; mtype = 4; end

%% Calculate
switch mtype
    case 1 % Average side length method
        edges = calc_schism_edge(Mobj);
        R = mean(edges, 2, 'omitnan');

    case 2 % Equal area circle radius method
        S = calc_schism_area(Mobj);
        R = sqrt(S/pi);

    case 3 % Square root area method
        S = calc_schism_area(Mobj);
        R = sqrt(S);

    case 4  % Circumcircle radius method
        edges = calc_schism_edge(Mobj);
        angles = calc_schism_angles(Mobj);
        inter_angle = angles(Mobj.i34==4, 2);

        L3 = edges(Mobj.i34==3, 1:3);
        L4 = edges(Mobj.i34==4, 1:4);
        Ld = sqrt(L4(:, 1).^2+L4(:, 2).^2-2*L4(:, 1).*L4(:, 2).*cosd(inter_angle));  % law of cosines
        L4_p1 = [L4(:, [1 2]) Ld(:)]; L4_p2 = [L4(:, [3 4]) Ld(:)];

        [~, S3, S4] = calc_schism_area(Mobj);

        R = nan(Mobj.nElems,1);
        R(Mobj.i34==3) = prod(L3, 2)./(4*S3);

        r4_p1 = prod(L4_p1, 2)./(4*S4(:,1));
        r4_p2 = prod(L4_p2, 2)./(4*S4(:,2));
        R(Mobj.i34==4) = max([r4_p1(:) r4_p2(:)], [], 2);  % use the largest circumscribed circle

    case 5 % Incircle radius method
        edges = calc_schism_edge(Mobj);
        angles = calc_schism_angles(Mobj);
        inter_angle = angles(Mobj.i34==4, 2);

        L3 = edges(Mobj.i34==3, 1:3);
        L4 = edges(Mobj.i34==4, 1:4);
        Ld = sqrt(L4(:, 1).^2+L4(:, 2).^2-2*L4(:, 1).*L4(:, 2).*cosd(inter_angle));  % law of cosines
        L4_p1 = [L4(:, [1 2]) Ld(:)]; L4_p2 = [L4(:, [3 4]) Ld(:)];

        [~, S3, S4] = calc_schism_area(Mobj);
        
        R = nan(Mobj.nElems,1);
        R(Mobj.i34==3) = 2*S3./sum(L3, 2);

        r4_p1 = 2*S4(:,1)./sum(L4_p1, 2);
        r4_p2 = 2*S4(:,2)./sum(L4_p2, 2);
        R(Mobj.i34==4) = r4_p1+r4_p2;   % using the sum of incircle radius
        
end
