function CR = calc_schism_cradius(Mobj)
% Calculate the circumradius (meters) on triangular mesh
%
%% Syntax
% CR = calc_schism_cradius(Mobj)
%
%% Description 
% CR = calc_schism_cradius(Mobj) calculates the circumradius (m) on the
% triangular mesh 
%
%% Example
% CR = calc_schism_cradius(Mobj);
%
%% Input Arguments
% Mobj --- the mesh object
%
%% Output Arguments
% CR --- circumradius of each triangular cell, can be regarded as the
% spatial resolution of mesh.
% 
%% Notes
% This function is adapted from the OceanMesh2D toolbox
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-05-18.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: calc_schism_area

%% Parse inputs
tri = Mobj.tri;
lon = Mobj.lon;
lat = Mobj.lat;

m_proj('Lambert', 'lon', [min(lon) max(lon)], 'lat', [min(lat) max(lat)])
[X,Y]= m_ll2xy(lon, lat);

% get the circumcenter radius of each element
TR = triangulation(tri, X,Y);
[~,cr] = circumcenter(TR);
% Get the element connectivity
[vtoe,nne] = VertToEle(tri);
% Make sure null value adds zero contribution
cr(end+1) = 0;
vtoe(vtoe == 0) = length(tri) + 1;
% Sum up all element contributions to each node and
% divide by number of connected elements
CR = sum(cr(vtoe))./nne;

% scale by earth radius
Re = 6378.137e3;
CR = Re*CR(:);

end

function [vtoe,nne] = VertToEle(t)
% this function is obtained from the OceanMesh2D toolbox

np = max(t(:)); ne = size(t,1);
nne = zeros(np,1);
for ie = +1 : ne                                                            %--Go through once and find max nnz.
    %for iv = +1 : +3
    %kjr loop unrolled.
    nm1=t(ie,1);
    nm2=t(ie,2);
    nm3=t(ie,3);
    %vrtx = t(ie,iv);
    nne(nm1,1) = nne(nm1,1) + 1;
    nne(nm2,1) = nne(nm2,1) + 1;
    nne(nm3,1) = nne(nm3,1) + 1;
    %end
end
mnz  = max(nne);                                                           %--max number of non-zeros
vtoe = zeros(mnz,np);                                                      %--vertex to element connectivity
nne = zeros(np,1);                                                         %--number of neighboring elements
for ie = +1 : ne
    nm1=t(ie,1);
    nm2=t(ie,2);
    nm3=t(ie,3);

    %kjr loop unrolled.
    nne(nm1,1) = nne(nm1,1) +1;
    nne(nm2,1) = nne(nm2,1) +1;
    nne(nm3,1) = nne(nm3,1) +1;

    vtoe(nne(nm1,1),nm1) = ie;
    vtoe(nne(nm2,1),nm2) = ie;
    vtoe(nne(nm3,1),nm3) = ie;

    %vtoe(nm1,nne(nm1,1)) = ie;
    %vtoe(nm2,nne(nm2,1)) = ie;
    %vtoe(nm3,nne(nm3,1)) = ie;
end
nne  = nne';
%vtoe = vtoe';

end