function hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag)
% Calculate the horizontal diffusivity
%
%% Syntax
% hdif = calc_schism_hdif(Mobj, gamma)
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max)
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag)
% 
%% Description
% hdif = calc_schism_hdif(Mobj, gamma) calculates the horizontal
%       diffusivity according to the mesh resolution.
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max) sepcifies the max. value.
% hdif = calc_schism_hdif(Mobj, gamma, hdif_max, disp_flag) determines
%       whether to display the results on the map.
% 
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store the mesh info.
% gamma - tuning factor; numeric
%       a tuning factor in estimating the horizontal diffusivity, ranging
%       in [0, 0.25]. Default: gamma = 0.15; 
% hdif_max - upper limit; numeric
%       the max. horizontal diffusivity. Default: hdif_max = 0.125;
% disp_flag --- display flag; char
%       the flat used to display the results (on/off). Default: disp_flag = 'on'. 
%
%% Output Arguments
% hdif - horizontal diffusivity; numeric
%       the horizontal diffusivity used to generate "hdif.gr3" file.
%
%% Notes
% This function aims to calculate the horizontal diffusivity for SCHISM
% model; In the formula, L_h is a length scale, equals to the minmum length
% of each triangle grid. gamma is a dimensionless constant between [0, 0.25)
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 26 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: calc_schism_bfric

%% Parse inputs
if nargin < 2; gamma = 0.15; end
if nargin < 3; hdif_max =  0.125; end
if nargin<4; disp_flag = 'on'; end
if gamma < 0 || gamma > 0.25; error('please make sure gamma is in [0, 0.25]'); end
if hdif_max>0.125; warning on; warning('hdif_max is greater than 1/8'); end

%% Calculate
edge_lens = calc_schism_edge(Mobj);
L_h = min(edge_lens, [], 2, 'omitnan');  % the shortest edge
L_h = (L_h-min(L_h))/(max(L_h)-min(L_h));
hdif_e = L_h.^2 * gamma/Mobj.dt;  % formula

hdif = convert_schism_var(Mobj, hdif_e, 'elem2node');
hdif = min(hdif, hdif_max);

disp(['hdif ranges from ', num2str(min(hdif)), ' to ', num2str(max(hdif))])
%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, hdif)
    hold on
    plot_schism_bnds(Mobj)
    colormap(jet(25))
    box on
    axis image
    auto_center
end
end
  
