function Mobj = gen_schism_SZ(Mobj, s_consts, zcors)
% Generate SZ (Sigma-Z) coordinates for SCHISM
% 
%% Syntax
% 
%
%% Description 
% 
%
%% Examples
%
%
%% Input Arguments
% Mobj --- mesh object
% zcors --- z-coordinate list, input as a vector; e.g. zcors = 20:5:100;
% the min. zcors is the transition depth between z-layers and sigma-layers.
% s_consts --- some consts for the sigma-formula. 
% s_consts = [h_c, theta_b, theta_f, ns]; h_c is a positive constant
% dictating the thickness of the bottom or surface layer that needs to be
% resolved, and theta_b and theta_f are constants that control the vertical
% resolution near the bottom and surface; ns is the # of sigma layers.
%
%% Output Arguments
% Mobj --- mesh object with 'depLayers' and 'vtype' loaded.
% 
%% Notes
% This function aims to generate the SZ (Sigma-Z) vertical grids for SCHISM
% model. The sigma-layers is on the top of z-layers. Purely-Z or -Sigma will
% be taken as a particular case.
%
% If the min. values in zcors > max. model depth, then a purely-S
% coordinate is adopted; if the min. values in zcors is zeros, a purely-Z
% is adopted
%
% If you want to create a set of 2D vertical grids, just use zcors = 100000,
% and s_consts = [100. 0. 1.e-4 2]; 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 
%
%% Parse inputs
if isfield(Mobj, 'vtype')
    error(['A vertical grid (',Mobj.vtype,') already exists!'])
end
%% Parse inputs
if nargin<2
    s_consts = [10, 0.7, 5, 20];
    zcors = 20:2:(ceil(max(Mobj.depth))+10);
end

if nargin<3
    zcors = 20:2:(ceil(max(Mobj.depth))+10);
end

h_c = s_consts(1);
theta_b = s_consts(2);
theta_f = s_consts(3);
ks = s_consts(4);
eta = 0;
%% SZ
zcors = sort(abs(zcors(:)));
if max(zcors)<max(Mobj.depth)
    error('Max. zcors must be greater than Max. model depth!')
end

h_s = min(zcors);
k = (1:ks)';
sigma = (k-1)/(1-ks);

Cs = @(sigma) (1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f) + theta_b*(tanh(theta_f*(sigma+0.5)) - tanh(theta_f/2))/(2*tanh(theta_f/2));

s_mas = nan(ks, Mobj.nNodes);

z_mas = repmat(zcors, 1, Mobj.nNodes);
for ii = 1:Mobj.nNodes
    h = Mobj.depth(ii);
    hw = min(h, h_s);
    hc0 = min(h_c, h);
    s_mas(:, ii) = eta*(1+sigma) + hc0*sigma + (hw-hc0)*Cs(sigma);
    ind = z_mas(:, ii)>h;
    ind_max = find(ind==0, 1, 'last' );
    z_mas(ind, ii) = nan;
    z_mas(ind_max+1, ii) = h;
end

ind_nans = sum(~isnan(z_mas), 2)==0;
z_mas(ind_nans,:) = [];

vtype = 'SZ';
sz_mas = [s_mas; -z_mas(2:end, :)];
if h_s > max(Mobj.depth)
    disp('Purely-S coordinate is adopted')
    sz_mas = s_mas;
    vtype = 'S';
    zcors = zcors(1);
end
if h_s == 0
    disp('Purely-Z coordinate is adopted')
    sz_mas = -z_mas;
    vtype = 'Z';
    sigma = [0; -1];
    s_consts = [100. 0. 1.e-4 2];
end

Mobj.vtype = vtype;
Mobj.s_conts = s_consts;
Mobj.sigma = sigma;
Mobj.zcors = -zcors(~ind_nans); % to avoid wrong # of layers when you have nans in z coordinates (Special thanks to Dr. Marcio Cintra for pointing this out)
Mobj.nLevs = sum(~isnan(sz_mas));
Mobj.maxLev = max(Mobj.nLevs);
Mobj.depLayers = sz_mas;

end


















