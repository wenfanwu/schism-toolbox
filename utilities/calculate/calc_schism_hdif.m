function hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag)
% calculate the horizontal diffusivity
%
%% Syntax
% hdif = calc_schism_hdif(Mobj, gamma)
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max)
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag)
%% Description
% hdif = calc_schism_hdif(Mobj, gamma) calculates the horizontal
% diffusivity according to the mesh resolution.
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max) sepcifies the max. value.
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag) determines
% whether to display the results on the map.
%% Input Arguments
% Mobj --- the mesh object
% gamma --- a tuning factor in estimating the horizontal diffusivity,
% ranging in [0, 0.25]. Default: gamma = 0.15;
% hdif_max --- max. horizontal diffusivity. Default: hdif_max = 6e-5;
% disp_flag --- flags to display the results (on/off). Default: disp_flag = 'off'. 
%
%% Output Arguments
% hdif --- horizontal diffusivity.
%
%% Notes
% this function aims to calculate the horizontal diffusivity for SCHISM
% model; In the formula, L_h is a length scale, equals to the minmum length
% of each triangle grid. gamma is a dimensionless constant between [0,0.25)
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2 Dec. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: disp_schism_var

%% Parse inputs
if nargin < 2
    gamma = 0.15;
end
if nargin < 3
    hdif_max = 6e-5;
end
if nargin<4
    disp_flag = 'off';
end

if gamma < 0 || gamma > 0.25
    error('please make sure gamma is in [0, 0.25]')
end
%% Calculate
L_h = zeros(1, Mobj.nElems);
for ii = 1:Mobj.nElems
    lonList = Mobj.lon(Mobj.tri(ii,:));
    latList = Mobj.lat(Mobj.tri(ii,:));
    L_h(ii) = min(hypot(lonList-lonList([2 3 1]), latList-latList([2 3 1])));
end
hdif_e = L_h.^2 * gamma/Mobj.dt;  % formula

nne = zeros(Mobj.nNodes,1);
indel = [];
for ii = 1:Mobj.nElems
    for jj = 1:3
        nd = Mobj.tri(ii, jj);
        nne(nd) = nne(nd)+1;
        indel(nne(nd), nd) = ii; %#ok<*AGROW>
    end
end

hdif = zeros(Mobj.nNodes,1);
for ii=1:Mobj.nNodes
    for jj =1:nne(ii)
        ie = indel(jj,ii);
        hdif(ii) = hdif(ii)+hdif_e(ie)/nne(ii);
    end
    hdif(ii) = min(hdif(ii), hdif_max);
end
hdif = min(hdif, hdif_max);
%% Display
disp(['hdif ranges from ', num2str(min(hdif)), ' to ', num2str(max(hdif))])
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, hdif)
    colormap(jet(25))
    box on
    axis image
end
end
  
