function elem_area = calc_schism_area(Mobj, ctype)
% Calculate the element area (m^2)
%
%% Syntax
% elem_area = calc_schism_area(Mobj)
% elem_area = calc_schism_area(Mobj, ctype)
% 
%% Description
% elem_area = calc_schism_area(Mobj) calculates the areas of elements on a
% triangulation mesh (units: m^3).
% elem_area = calc_schism_area(Mobj, ctype) determines the coorinate type
% (geographic/cartesian).
%
%% Example
% elem_area = calc_schism_area(Mobj, 'geographic');
% 
%% Input Arguments
% Mobj --- the mesh object
% ctype --- coordinate type (geographic/cartersian); Default: ctype = 'geographic'; 
%
%% Output Arguments
% elem_area --- the areas of elements  (units: m^2).
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 8 Mar. 2022. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: distance

%% Parse inputs
if nargin < 2
    ctype = 'geographic';
end
lon_verts = Mobj.lon(Mobj.tri);
lat_verts = Mobj.lat(Mobj.tri);
nElems = size(Mobj.tri,1);

%% Calculation
if strncmpi(ctype, 'geographic', 3)
    disp('calculate the element areas on a geographic coordinate')
    N = pi/180*6378000;
    a = distance(lat_verts(:, 1), lon_verts(:, 1), lat_verts(:, 2), lon_verts(:, 2))*N;
    b = distance(lat_verts(:, 2), lon_verts(:, 2), lat_verts(:, 3), lon_verts(:, 3))*N;
    c = distance(lat_verts(:, 1), lon_verts(:, 1), lat_verts(:, 3), lon_verts(:, 3))*N;
    p = sum([a(:), b(:), c(:)], 2)/2;

    elem_area = sqrt(abs(p.*(a-p).*(p-b).*(p-c))); % Haron Formula
else
    disp('calculate the element areas on a Cartesian coordinate')
    side_vec1 = [zeros(nElems,1) lon_verts(:,2)-lon_verts(:,1) lon_verts(:,3)-lon_verts(:,1)];
    side_vec2 = [zeros(nElems,1) lat_verts(:,2)-lat_verts(:,1) lat_verts(:,3)-lat_verts(:,1)];
    
    elem_area = sum(cross(side_vec1, side_vec2)/2, 2);  % Cross-product
end

end