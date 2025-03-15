function Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, s_consts, dz_bot_min)
% Generate LSC2-type vertical grids
% 
%% Syntax
% Mobj = gen_schism_LSC2(Mobj)
% Mobj = gen_schism_LSC2(Mobj, dep_edges)
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers)
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers, s_consts)
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers, s_consts, dz_bot_min)
%
%% Description 
% Mobj = gen_schism_LSC2(Mobj) generates LSC2-stype vertical layers
% Mobj = gen_schism_LSC2(Mobj, dep_edges) specifies the depth edges
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers) specifies the layers of each depth bin
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers, s_consts) specifies the streching constants
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_layers, s_consts, dz_bot_min) specifies the minimum BBL thickness.
% 
%% Example
% dep_edges = [10, 20, 30, 45, 55, 65, 75, 90];
% dep_nums =  [20 21 22 23 24 25 27 28];
% Mobj = gen_schism_LSC2(Mobj, dep_edges, dep_nums, [4 5 3 5], 0.25);
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       a datastruct storing the grid info.
% dep_edges - depth edges; double
%       depth edges used to partition the depth bins.
% dep_nums - the # of layers; double
%       the # of depth layers for each depth bin.
% s_const - stretching constants. double
%       the stretching constants used for sigma-function
%       s_const = [VSTRETCHING, rtheta_s, rtheta_b, TCLINE]; 
%       VSTRETCHING determines the stretching function used here.
%       rtheta_s/rtheta_b determines the weight of the refinement in
%       surface/bottom layer. s_const can tweak the sparsity of vertical
%       grids, so as to realize the functionality of surface/bottom layer
%       refinement; default: s_const = [4 5 3 5];
%       here are some examples for s_const: 
%       surface refinement                     bottom refinement 
%       s_consts = [4 5 3 5];                   s_consts = [3 1 3 5]; 
%       s_consts = [2 7 0.1 5];          
%       s_consts = [2 9 0.1 5];          
% dz_bot_min - minimum thickness; double
%       the minmum thickness (meters) of the bottom bounday layer; 
%       default: dz_bot_min = 0.35.
%
%% Output Arguments
% Mobj - the mesh object; datastruct
%       the mesh object with vertical layers added.
%
%% Notes
% This function aims to generate the LSC2 vertical grids for SCHISM
% model. It is adapted from the "gen_vqs_Rutgers.f90" in the SCHISM
% source code. Refer to the original FORTRAN code for more details.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last updated on 25 Feb 2025. 
% Email: wwu@vims.edu
% 
% See also: gen_schism_SZ and write_schism_vgrid

%% Parse inputs
if isfield(Mobj, 'vtype')
    error(['vertical grids (',Mobj.vtype,') already existed!'])
end

global TCLINE VSTRETCHING rtheta_s rtheta_b np ne xnd ynd dp hsm m_vqs nv_vqs %#ok<*GVMIS> 

if nargin < 2
    dep_edges = unique(ceil(prctile(Mobj.depth, [20 30 40 50 60 70 80 85 90 95 98 100])));
end
if isempty(dep_edges)
    dep_edges = unique(ceil(prctile(Mobj.depth, [20 30 40 50 60 70 80 85 90 95 98 100])));
end
m_vqs= numel(dep_edges);

if nargin < 3
    dep_nums(1:m_vqs) = 18+(1:m_vqs);
end
if isempty(dep_nums)
    dep_nums(1:m_vqs) = 18+(1:m_vqs);
end

if nargin < 4
    s_consts = [4 5 3 5];
end
if nargin <5
    dz_bot_min = 0.35;
end

hsm = dep_edges;
nv_vqs = dep_nums;

%% ================================================
%        ================ Streching const. =============
% ================================================
% Notes: rtheta_s and rtheta_b are the surface and bottom control
% parameters.  0 < rtheta_s ≤ 20 and 0 ≤ rtheta_b ≤ 1.
% 1) The larger the value of rtheta_s, the more resolution is kept above h_c.
% 2) For rtheta_b = 0, the resolution all goes to the surface as rtheta_s is increased.
% 3) For rtheta_b = 1, the resolution goes to both the surface and the
% bottom equally as rtheta_s is increased. 
% 4) For rtheta_s ≠ 0, there is a subtle mismatch in the discretization of the
% model equations, for instance in the horizontal viscosity term. We
% recommend that you stick with reasonable values of rtheta_s, say rtheta_s ≤ 8.  
% 5) Some applications turn out to be sensitive to the value of rtheta_s used.

VSTRETCHING = s_consts(1); % different stretching function
rtheta_s = s_consts(2);
rtheta_b = s_consts(3);
TCLINE  = s_consts(4);
% ================================================
% ================================================
% ================================================

ne = Mobj.nElems;
np = Mobj.nNodes;
xnd = Mobj.lon;
ynd = Mobj.lat;
dp = Mobj.depth;
nvrt_m = nv_vqs(m_vqs);
nvrt=nvrt_m;

sigma_vqs = nan(nvrt,np);

if m_vqs<2
    error('Check vgrid.in: m_vqs')
end
if hsm(1)<0
    error('hsm(1)<0')
end
for m=2:m_vqs
    if hsm(m)<=hsm(m-1)
        error(['Check hsm: ', num2str(m), num2str(hsm(m)), num2str(hsm(m-1))])
    end
end
%    Other consts.
%    Stretching const. for the 1st master grid and also for depth <= hsm(1)
%     |a_vqs0|<=1 (1: skew toward bottom; -1: toward surface; 0: no bias)
a_vqs0 = -1;
%     Generate a master vgrid (z_mas)
etal=0;    %used in master grid only; elev.
if etal<=-hsm(1)
    error(['elev<hsm', num2str(etal)])
end

disp(['nvrt in master vgrid=', num2str(nvrt_m)])
z_mas = -1.e5*ones(nvrt_m, m_vqs);

opt_val = 1;
a_vqs =zeros(m_vqs,1);
for m=1:m_vqs
    for k=1:nv_vqs(m)
        sigma=(k-1.0)/(1.0-nv_vqs(m));
        % Alternative transformations below
        % Option 1: quadratic
        if opt_val==0
            a_vqs(m)=max(-1.d0,a_vqs0-(m-1)*0.03);
            tmp=a_vqs(m)*sigma*sigma+(1+a_vqs(m))*sigma; %transformed sigma
            z_mas(k,m)=tmp*(etal+hsm(m))+etal;
        end
        % Option 2: S
        if opt_val==1
            theta_b=0.0;
            theta_f=4.0; % 0 - 20 , toward 0 like a traditionalsigma with uniform spacing
            
            cs1 = (1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f);
            cs2 = theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5);
            cs=cs1+cs2;
            z_mas(k,m)=etal*(1+sigma)+hsm(1)*sigma+(hsm(m)-hsm(1))*cs;
        end
    end
    % Option 3: Rutgers Coordinate : z_mas above will be erased
    if opt_val==1
        k = nv_vqs(m);
%         sc_w = 0.0; Cs_w = 0.0;
        [Cs_w, sc_w] = SIGMA_RUTGERS(k);
        %            !Compute the sigma coordinate from 0 to hsm(m)
        V_RUTGERS = SIGMA_RUTGERS_VEC(hsm(m),k,sc_w,Cs_w);
        %             !WRITE(*,*) V_RUTGERS
        z_mas(1:k,m) = V_RUTGERS(1:k)*hsm(m);
    end
end

% Read in hgrid
i34 = 3*ones(ne,1);
elnode = Mobj.tri';
eta2=etal*ones(np,1);

% Check max depth
dpmax=max(dp);
if dpmax>hsm(m_vqs)
    error(['Max depth exceeds master depth:', num2str(dpmax), num2str(hsm(m_vqs))])
end

ns = Mobj.nEdges;
isidenode = Mobj.edg';
disp(['ns=', num2str(ns)])

% Preparation
% --------------
%  Shallow area with depth<= hsm
%  Get Stretching Coeff sc_w and Cs_w
[Cs_w,sc_w] = SIGMA_RUTGERS(nv_vqs(1));

% Compute the sigma coordinate on the whole grid a la RUTGERS,
% from 0 to hsm : S_RUTGERS is Updated to be used later line 303
% in computing "sigma_vqs" and "znd"
S_RUTGERS = SIGMA_RUTGERS_MAT(nv_vqs(1), sc_w,Cs_w);

% ---------------------------------------------------
% Compute zcoor
% ---------------------------------------------------
znd=-1.e6*ones(nvrt,np);
m0 = zeros(np,1);
kbp = zeros(np,1);
bbb=1;
for i=1:np
    % ----------------------------------------------------------
    % Shallow region
    % -----------------------------------------------------------
    if dp(i)<=hsm(1) %shallow
        kbp(i)=nv_vqs(1);
        for k=1:nv_vqs(1)
            if bbb==0 % Original Way, The Shism way
                sigma=(k-1.0)/(1.0-nv_vqs(1));
                sigma_vqs(k,i)=a_vqs0*sigma*sigma+(1+a_vqs0)*sigma; %transformed sigma
                znd(k,i)=sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i);
            else  %New way to compute sigma in shallow using Rutgers Streching Function
                sigma_vqs(k,i)=S_RUTGERS(k,i);
                znd(k,i) = sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i);
            end
        end
        continue
    end
    % ----------------------------------------------------------
    % Deep region : If RUTGERS TRUE Line 188, then z_mas depths are those
    %                       computed with Rutgers Stretching function
    % ----------------------------------------------------------
    %Find a master vgrid
    for m=2:m_vqs
        if dp(i)>hsm(m-1) && dp(i)<=hsm(m)
            m0(i)=m;
            zrat=(dp(i)-hsm(m-1))/(hsm(m)-hsm(m-1)); %(0,1]
            break
        end
    end %m
    if m0(i)==0
        error(['Failed to find a master vgrid:', num2str(i) num2str(dp(i))])
    end

    kbp(i) = 0;
    for k=1:nv_vqs(m0(i))
        z1=z_mas(min(k,nv_vqs(m0(i)-1)),m0(i)-1);
        z2=z_mas(k,m0(i));
        z3=z1+(z2-z1)*zrat;
        if z3>=-dp(i)+dz_bot_min
            znd(k,i)=z3;
        else
            kbp(i)=k;
            break
        end
    end %k
    
    if kbp(i)==0
        error(['Failed to find a bottom:',num2str(i),' ',num2str(dp(i)),' ',num2str(z3), ' ',num2str(z_mas(1:nv_vqs(m0(i)))), ' ',num2str(m0(i))])
    end
    znd(kbp(i),i)=-dp(i);
    %Check order
    for k=2:kbp(i)
        if znd(k-1,i)<=znd(k,i)
            error(['Inverted z:',num2str(i), ' ',num2str(dp(i)), ' ',num2str(m0(i)), ' ',num2str(k), ' ',num2str(znd(k-1,i)), ' ',num2str(znd(k,i))])
        end
    end %k
end %i=1,np

%    Extend beyond bottom for plotting
for i=1:np
    znd(kbp(i)+1:nvrt,i)=-dp(i);
end %i

if 1==2
    % !----------------------------
    % !     Check case where a side saddles hsm(1) - try to make the # of level equal
    % !Error: need to iterate?
    for i=1:ns
        n1=isidenode(1,i);
        n2=isidenode(2,i);
        while kbp(n1)==kbp(n2) || (dp(n1)-hsm(1))*(dp(n2)-hsm(1))>0
            disp(['Correcting layers near hsm(1)...', num2str(n1),num2str(n2),num2str(real(xnd(n1))),num2str(real(xnd(n2)))])
        end
        if kbp(n1)>kbp(n2)
            in1=n2; in2=n1;
        else
            in1=n1; in2=n2;
        end
        
        k_add=kbp(in2)-kbp(in1);
        if k_add<=0
            error('k_add<=0')
        end
        for k=1:k_add
            znd(kbp(in1)-1+k,in1)=znd(kbp(in1)-1,in1)-k*(znd(kbp(in1)-1,in1)+dp(in1))/(k_add+1);
        end %k
        znd(kbp(in1)+k_add,in1)=-dp(in1);
        kbp(in1)=kbp(in1)+k_add;
        
        %Check
        for k=2:kbp(in1)
            if znd(k-1,in1)<=znd(k,in1)
                error(['Inverted z (2):',num2str(n1),num2str(n2),num2str(dp(n1)),num2str(dp(n2)), ...
                    num2str(k),num2str(znd(k-1,in1)),num2str(znd(k,in1))])
            end
        end %k
    end %i=1,ns
end %1==2

if min(kbp)<2
    disp(['# of levels <2:', min(kbp)])
    error('Check parameters like dz_bot_min etc')
end

nvrt=max(kbp);
disp(['Final nvrt=', num2str(nvrt)])
%     # of prisms
nprism=0;
for i=1:ne
    kbpl=max(kbp(elnode(1:i34(i),i)));
    nprism=nprism+kbpl;
end %i
disp(['# of prisms=', num2str(nprism)])
% disp(['Average # of layers', num2str(nprism/ne)])

% Output in SCHISM convention
% fid = fopen([Mobj.aimpath, 'vgrid.in'],'wt');
% fprintf(fid, '%d !ivcor\n', 1);
% fprintf(fid, '%d\n', nvrt);
for i=1:np
    if dp(i)<=hsm(1)
        %  sigma_vqs already assigned
    else
        sigma_vqs(1,i)=0;
        sigma_vqs(kbp(i),i)=-1;
        for k=2:(kbp(i)-1)
            sigma_vqs(k,i)=(znd(k,i)-eta2(i))/(eta2(i)+dp(i));
        end %k
    end
%  Check order
    for k=2:kbp(i)
        if sigma_vqs(k,i)>=sigma_vqs(k-1,i)
            error(['Inverted sigma:',num2str(i), num2str(k), num2str(dp(i)), num2str(sigma_vqs(k,i)), num2str(sigma_vqs(k-1,i))])
        end
    end
%     nLevs = nvrt+1-kbp(i);
%     sigLevs = rot90(sigma_vqs(1:kbp(i),i),3);
%     formatStr = ['%d %8d ', repmat('% 14.6f', 1, length(sigLevs)), '\n'];
%     fprintf(fid, formatStr, i, nLevs, sigLevs);
end %i
% fclose(fid);

Mobj.vtype = 'LSC2';
Mobj.nLevs = sum(~isnan(sigma_vqs));
Mobj.maxLev = max(Mobj.nLevs);
disp(['The mean # of levels is ', num2str(mean(Mobj.nLevs), '%.2f')])

Mobj.master_grid = 'unknown';
vgrids = sigma_vqs(1:Mobj.maxLev,:);

Mobj.depLayers = bsxfun(@times, Mobj.depth', vgrids);
Mobj.vgrids = vgrids;

end

function [Cs_w, sc_w] = SIGMA_RUTGERS(KB)  % S_RUTGERS
% !
% ! JEROME (IRD Noumea, 5-Oct-2017)
% !
% ! Support to use several Vertical Stretching function for vertical
% ! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
% ! ------------------------------------------------------------------
% ! Various possible vertical stretching are provided. They are tunable
% ! by using different value for Vstretching (2, 3 or 4), rtheta_s, rtheta_b
% ! and Hc (Tcline). The original vertical stretching function from Song
% ! and
% ! Haidvogel (1994) is not supplied, but can be approached by setting
% ! Vstretching = 2
% !
% ! See: Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic
% ! modeling system (ROMS): a split-explicit, free-surface,
% ! topography-following-coordinate oceanic model, Ocean
% ! See details: https://www.myroms.org/wiki/Vertical_S-coordinate
% !
% ! In Subroutine SIGMA_RUTGERS :
% ! the original ROMS/RUTGERS/UCLA algorithm to compute sigma coordinate S
% ! The vertical stretching function Cs at ROMS (SHISM) W-points (layer
% ! interfaces) are applied.
% !
% ! In Subroutine SIGMA_RUTGERS_MAT
% ! Vertical depths through the whole Vertices and levels are applied
% ! using the Rutgers Vertical streching function.
% ! In Section 3, a wrapper is applied to set the value of sigma
% ! coordinates compatible with FVCOM/SCHISM
%
% ! In Subroutine SIGMA_RUTGERS_VEC
% ! Vertical depths at discrete vertices (in deep region) are applied
% ! using the Rutgers Vertical streching function.
% ! In Section 3, a wrapper is applied to set the value of sigma
% ! coordinates compatible with FVCOM/SCHISM

% ! ------------------------------------------------------------------
% ! Section 1 : Sigma Vertical coordinate Definition by using various
% ! stretching function as introduced by the RUTGERS/UCLA team
% ! --------------------------------------------------------------------
% ! Set S(sigma(k)), the non dimensionnal stretched vertical coordinate
% !  with :  -1 <= S <= 0  ;  S= 0 at the surface   S = -1 at the bottom
% ! Set  C : nondimensional vertical stretching function, C(sigma(k)),
% !  with :  -1 <= C(s) <= 0  : C= 0 at the surface   C = -1 at the bottom
% !
% ! Local naming convention (ROMS Rutgers Like) :
% ! sc_w at W-point i.e layer interface (dimension KB)
% ! sc_r at Rho-point i.e mid layer     (dimension KBM1)
% ! Cs_w at W-point i.e layer interface (dimension KB)
% ! Cs_r at Rho-point i.e mid layer     (dimension KBM1)
% ! --------------------------------------------------------------------
global VSTRETCHING rtheta_s rtheta_b

% Wrapper :
KBm1 = KB-1;
% !-----------------------------------------------------------------------
% ! Vstretching = 2 : The New  A. Shchepetkin new vertical stretching
% !                   function.
% ! See :
% !    Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic !
% !         modeling system (ROMS): a split-explicit, free-surface,      !
% !         topography-following-coordinate oceanic model, Ocean         !
% !         Modelling, 9, 347-404.
% !-----------------------------------------------------------------------
sc_w = zeros(KB,1);
Cs_w = zeros(KB,1);

if VSTRETCHING==2
    Aweight = 1.0;
    Bweight = 1.0;
    ds=1.0/KBm1;
    
    % ! W-Point layer Interface
    sc_w(KB)=0.0;
    Cs_w(KB)=0.0;
    for k=(KBm1-1):-1:1
        cff_w  = ds*(k-KBm1);
        sc_w(k+1)= cff_w;
        if rtheta_s>0
            Csur=(1.0-cosh(rtheta_s*cff_w))/(cosh(rtheta_s)-1.0);
            if rtheta_b>0
                Cbot=sinh(rtheta_b*(cff_w+1.0))/sinh(rtheta_b)-1.0;
                Cweight=(cff_w+1.0)^Aweight*(1.0+(Aweight/Bweight)*(1.0-(cff_w+1.0)^Bweight));
                Cs_w(k+1)=Cweight*Csur+(1.0-Cweight)*Cbot;
            else
                Cs_w(k+1)=Csur;
            end
        else
            Cs_w(k+1)=cff_w;
        end
    end
    sc_w(1) = -1;
    Cs_w(1) = -1;
    % !-----------------------------------------------------------------------
    % ! Vstretching = 3 : R. Geyer stretching function for high bottom
    % !                   boundary layer resolution
    % !-----------------------------------------------------------------------
elseif VSTRETCHING==3
    exp_sur=rtheta_s;
    exp_bot=rtheta_b;
    Hscale=3.0;
    ds=1.0/(KBm1);
    % ! W-Point layer Interface
    sc_w(KB)=0.0;
    Cs_w(KB)=0.0;
    for k=(KBm1-1):-1:1
        cff_w  = ds*(k-KBm1);
        sc_w(k+1)= cff_w;
        Cbot= log(cosh(Hscale*(cff_w+1.0)^exp_bot))/log(cosh(Hscale))-1.0;
        Csur=-log(cosh(Hscale*abs(cff_w)^exp_sur))/log(cosh(Hscale));
        Cweight=0.5*(1.0-tanh(Hscale*(cff_w+0.5)));
        Cs_w(k+1)=Cweight*Cbot+(1.0-Cweight)*Csur;
    end
        sc_w(1) = -1;
    Cs_w(1) = -1;
    % -----------------------------------------------------------------------
    %  Vstretching = 4 : A. Shchepetkin (UCLA-ROMS, 2010) double vertical
    %                    stretching function
    % -----------------------------------------------------------------------
elseif VSTRETCHING==4
    ds=1.0/KBm1;
    % ! W-Point layer Interface
    sc_w(KB)=0.0;
    Cs_w(KB)=0.0;
    
    for k=(KBm1-1):-1:1
        cff_w = ds*(k-KBm1);
        sc_w(k+1)  = cff_w;
        if rtheta_s>0.0
            Csur=(1.0-cosh(rtheta_s*cff_w))/(cosh(rtheta_s)-1.0);
        else
            Csur = (cff_w^2)*(-1.0);
        end
        if rtheta_b>0.0
            Cbot=(exp(rtheta_b*Csur)-1.0)/(1.0-exp(-rtheta_b));
            Cs_w(k+1)=Cbot;
        else
            Cs_w(k+1)=Csur;
        end
    end
    sc_w(1) = -1;
    Cs_w(1) = -1;
    
    % ! ----------------------
else
    error('SIGMA_RUTGERS: Wrong value for VSTRETCHING, only 2,3 or 4 allowed')
end
end

%--------------------------------------------------------------------

function S_RUTGERS = SIGMA_RUTGERS_MAT(KB,sc_w,Cs_w)
% !
% ! JEROME (IRD Noumea, 5-Oct-2017)
% !
% ! Support to use several schemes for generalized vertical
% ! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
% ! ------------------------------------------------------------------

% ! ------------------------------------------------------------------
% ! Section 2 : Compute Vertical Height as in ROMS RUTGERS/UCLA
% !-----------------------------------------------------------------------
% !  New formulation: Compute vertical depths (meters, negative) at
% !                   RHO- and W-points, and vertical grid thicknesses.
% !  Various stretching functions are possible, as defined above.
% !
% !         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
% !
% !         Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
% !
% !         but with zeta = 0
% !-----------------------------------------------------------------------
global TCLINE np dp hsm


KBm1 = KB-1;
hc = TCLINE;

z_w = zeros(np,KB);
H = 0.0;

for  I=1:np
    H(I) = min(dp(I),hsm(1));
    z_w(I,1) = -1.0;
end

for k=1:KBm1
    cff_w  = hc*sc_w(k+1);
    cff1_w = Cs_w(k+1);
    for I=1:np
        hwater=H(I);
        hinv=1.0/(hc+hwater);
        cff2_w=(cff_w+cff1_w*hwater)*hinv;
        z_w(I,k+1)=cff2_w;
    end
end

%  ------------------------------------------------------------------
%  Section 3 : WRAPPER : ROMS vert. coord. to SCHISM vert. coord.
%  -----------------------------------------------------------------------
S_RUTGERS = zeros(KB,np);
for I=1:np
    for K=1:KB
        KK=KB-K+1;
        S_RUTGERS(K,I) = z_w(I,KK);
    end
end
end

function V_RUTGERS = SIGMA_RUTGERS_VEC(H,KB,sc_w,Cs_w)
% ! JEROME (IRD Noumea, 5-Oct-2017)
% !
% ! Support to use several schemes for generalized vertical
% ! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
% ! ----------------------------------------------------------------
% ! ------------------------------------------------------------------
% ! Section 2 : Compute Vertical Height as in ROMS RUTGERS/UCLA
% !-----------------------------------------------------------------------
% !  New formulation: Compute vertical depths (meters, negative) at
% !                   RHO- and W-points, and vertical grid thicknesses.
% !  Various stretching functions are possible, as defined above.
% !
% !         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
% !
% !         Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
% !
% !         but with zeta = 0
% !-----------------------------------------------------------------------

global TCLINE

hc = TCLINE;
KBm1 = KB-1;

z_w = zeros(KB,1);
z_w(1) = -1.0;
for k=1:KBm1
    cff_w = hc*sc_w(k+1);
    cff1_w = Cs_w(k+1);
    
    hwater=H;
    hinv=1.0/(hc+hwater);
    cff2_w=(cff_w+cff1_w*hwater)*hinv;
    z_w(k+1)=cff2_w;
end
% ! ------------------------------------------------------------------
% ! Section 3 : WRAPPER : ROMS vert. coord. to SCHISM vert. coord.
% !-----------------------------------------------------------------------
V_RUTGERS = zeros(KB,1);
for K=1:KB
    KK=KB-K+1;
    V_RUTGERS(K) = z_w(KK);
end
end


