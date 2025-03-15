function Mobj = gen_schism_SZ(Mobj, s_consts, zcors)
% Generate SZ-type (Sigma-Z) vertical grids
% 
%% Syntax
% Mobj = gen_schism_SZ(Mobj)
% Mobj = gen_schism_SZ(Mobj, s_consts)
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors)
%
%% Description 
% Mobj = gen_schism_SZ(Mobj) generates sigma-z vertical grids
% Mobj = gen_schism_SZ(Mobj, s_consts) specifies the streching constants
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors) specifies the z-coordinates
%
%% Examples
% s_consts = [10, 0.7, 5, 20];
% zcors = 20:2:(fix(max(Mobj.depth))+10);
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors); % sigma-z hybrid
%
% s_consts = [10, 0.7, 5, 20];
% zcors = max(Mobj.depth)+100;
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors); % purely-sigma
%
% s_consts = [10, 0.7, 5, 20]; 
% zcors = 0:2:(ceil(max(Mobj.depth))+10);
% Mobj = gen_schism_SZ(Mobj, s_consts, zcors); % purely-z 
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       a datastruct storing the grid info.
% s_consts - streching constants; double
%       some consts for the sigma-formula. 
%       s_consts = [h_c, theta_b, theta_f, ns]; h_c is a positive constant
%       dictating the thickness of the bottom or surface layer that needs
%       to be resolved, and theta_b and theta_f are constants that control
%       the vertical resolution near the bottom and surface; ns is the # of
%       sigma layers. default; s_consts = [10, 0.7, 5, 20];
% zcors - z-coordinates; double
%       the vertical coordinates used for z-layers; e.g. zcors = 20:5:100;
%       the min. zcors is the transition depth between z-layers and
%       sigma-layers. default: zcors = 20:2:(ceil(max(Mobj.depth))+10);
%
%% Output Arguments
% Mobj - the mesh object; datastruct
%       the mesh object with vertical layers added.
% 
%% Notes
% As for sigma-z vertical grids, the sigma-layer is on the top of z-layers,
% and purely-z or -sigma grids will be taken as particular cases.
%
% If the min. values in zcors > max. model depth, then a purely-s
% grid is used; if the min. value in zcors is zero, a purely-z grid will be
% used instread.
%
% If you want to create a set of 2D vertical grids, just use zcors = 100000,
% and s_consts = [100. 0. 1.e-4 2]; 
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last updated on 25 Feb 2025. 
% Email: wwu@vims.edu
% 
% See also: gen_schism_LSC2 and write_schism_vgrid

%% Parse inputs
if isfield(Mobj, 'vtype')
    error(['vertical grids (',Mobj.vtype,') already existed!'])
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
    disp('purely-sigma coordinate is used')
    sz_mas = s_mas;
    vtype = 'S';
    zcors = zcors(1);
end
if h_s == 0
    disp('purely-z coordinate is used')
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


















