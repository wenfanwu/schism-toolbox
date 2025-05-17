function shapiro_val = calc_schism_shapiro(Mobj, slope_list, val_max, disp_flag)
% Calculate the shapiro filter based on bathymetric slope.
%
%% Syntax
% shapiro_val = calc_schism_shapiro(Mobj, slope_list)
% shapiro_val = calc_schism_shapiro(Mobj, slope_list, val_max)
% shapiro_val = calc_schism_shapiro(Mobj, slope_list, val_max, disp_flag)
% 
%% Description
% shapiro_val = calc_schism_shapiro(Mobj, slope_list)calculates the slope
%       filter for shapiro.gr3 file.
% shapiro_val = calc_schism_shapiro(Mobj, slope_list, val_max) specifies
%       the maximum filter strength.
% shapiro_val = calc_schism_shapiro(Mobj, slope_list, val_max, disp_flag)
%       specifies the display flags (on/off).
% 
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct used to store the mesh info.
% slope_list - threshold value for slope; numeric
%       a two-element vector for turning the slope filter. Default: slope_list = [0.001, 0.05];
% val_max - maximum shapiro values; numeric
%       the maximum filter strength. Default: val_max = 0.5;
%
%% Output Arguments
% shapiro_val - shapiro values; numeric
%       the shapiro values to be written.
%
%% Notes
% This function was directly translated from the function named
% gen_slope_filter2.f90 in SCHISM source code. Refers to the original
% FORTRAN script for more details.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 2 Dec. 2021. 
% Email: wwu@vims.edu
% 
% See also: write_schism_gr3

%% Parse inputs
if strncmpi(Mobj.coord, 'geographic', 3); ics = 2; else; ics = 1; end
if nargin < 2; slope_list = [0.001, 0.05]; end
if nargin < 3; val_max = 0.5; end
if nargin < 4; disp_flag = 'off'; end
if val_max > 0.5; error('the val_max can not be greater than 0.5 or less than 0'); end

%% Calculate
slope_min = slope_list(1);
threshold_slope = slope_list(2);
ne = Mobj.nElems;
np = Mobj.nNodes;

mnei=25;
rearth_eq=6378206.4;
rearth_pole=6378206.4;
nwild = zeros(3,1);
area = zeros(np,1);
dldxy = zeros(2,3);
slope = zeros(ne,1);
nne = zeros(np,1);
indel = zeros(mnei,np);
shapiro_val = zeros(np,1);
% hdif_e = zeros(ne,1);
% rlh = zeros(4,1);
lframe = zeros(3,3);
xnd = zeros(np,1);
ynd = zeros(np,1);
znd = zeros(np,1);
xloc = zeros(3,1);
yloc = zeros(3,1);

x = Mobj.lon;
y = Mobj.lat;
dp = Mobj.depth;

for i=1:np
    if ics==2
        xtmp=x(i)/180*pi;
        ytmp=y(i)/180*pi;
        xnd(i)=rearth_eq*cos(ytmp)*cos(xtmp);
        ynd(i)=rearth_eq*cos(ytmp)*sin(xtmp);
        znd(i)=rearth_pole*sin(ytmp);
    end
end

i34 = Mobj.i34;
elnode = Mobj.tri';

for i=1:ne
    %skip shallow elem
    hmax=max(dp(elnode(1:i34(i),i)));
    if hmax<=100
        %         cycle
    end
    for m=1:(i34(i)-2) % split into tri
        if m==1
            nwild(1:3)=[1 2 3];
        else %quad
            nwild(1:3)=[1 3 4];
        end %m
        
        n1=elnode(nwild(1),i);
        n2=elnode(nwild(2),i);
        n3=elnode(nwild(3),i);
        
        signa = @(x1,x2,x3,y1,y2,y3) ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2;
        if ics==1
            area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3));
        else %lon/lat
            %Construct local frame with 1,2 as local x-axis, and 2,1)x(3,1) as local z.
            %Compute local z-axis 1st. [lframe(1:3,1:3): 1st index is
            %component]
            lframe(1,1)=xnd(n2)-xnd(n1);
            lframe(2,1)=ynd(n2)-ynd(n1);
            lframe(3,1)=znd(n2)-znd(n1);
            
            [lframe(1,3),lframe(2,3),lframe(3,3)] = cross_product(lframe(1,1),lframe(2,1),lframe(3,1), xnd(n3)-xnd(n1),ynd(n3)-ynd(n1),znd(n3)-znd(n1));
            [lframe(1,2),lframe(2,2),lframe(3,2)] = cross_product(lframe(1,3),lframe(2,3),lframe(3,3), lframe(1,1),lframe(2,1),lframe(3,1));
            
            xloc(1)=0; yloc(1:2)=0;
            for j=1:3
                rnorm=sqrt(lframe(1,j)^2+lframe(2,j)^2+lframe(3,j)^2);
                if rnorm==0
                    error('0 vector')
                end
                lframe(:,j)=lframe(:,j)/rnorm;
                if j==1
                    xloc(2)=rnorm;
                end
            end
            
            xloc(3)=(xnd(n3)-xnd(n1))*lframe(1,1)+(ynd(n3)-ynd(n1))*lframe(2,1)+ (znd(n3)-znd(n1))*lframe(3,1);
            yloc(3)=(xnd(n3)-xnd(n1))*lframe(1,2)+(ynd(n3)-ynd(n1))*lframe(2,2)+ (znd(n3)-znd(n1))*lframe(3,2);
            
            area(i)=signa(xloc(1),xloc(2),xloc(3),yloc(1),yloc(2),yloc(3));
        end %ics
        
        if area(i)<=0
            error(['Negative area at', num2str(i)])
        end
        
        for j=1:3
            nj1=j+1;
            nj2=j+2;
            if nj1>3
                nj1=nj1-3;
            end
            if nj2>3
                nj2=nj2-3;
                nd1=elnode(nwild(nj1),i);
                nd2=elnode(nwild(nj2),i);
            end
            if ics==1
                dldxy(1,j)=(y(nd1)-y(nd2))/2/area(i);
                dldxy(2,j)=(x(nd2)-x(nd1))/2/area(i);
            else
                dldxy(1,j)=(yloc(nj1)-yloc(nj2))/2/area(i);
                dldxy(2,j)=(xloc(nj2)-xloc(nj1))/2/area(i);
            end %ics
        end%j
        slx=sum(dp(elnode(nwild(1:3),i)).*dldxy(1,:)', 'omitnan');
        sly=sum(dp(elnode(nwild(1:3),i)).*dldxy(2,:)', 'omitnan');
        tmp=sqrt(slx.^2+sly.^2);
        if tmp>slope_min
            slope(i)=max(slope(i),tmp);
        end
    end %m
end %i=1,ne

% Neighborhood
for i=1:ne
    for j=1:i34(i)
        nd=elnode(j,i);
        nne(nd)=nne(nd)+1;
        if nne(nd)>mnei
            error(['Too many neighbors', num2str(nd)])
        end
        indel(nne(nd),nd)=i;
    end
end

for i=1:np
    slopemax=0;
    for j=1:nne(i)
        ie=indel(j,i);
        slopemax=max(slopemax, slope(ie));
    end %j
    shapiro_val(i)=val_max*tanh(2*slopemax/threshold_slope);
    if shapiro_val(i)<=1.e-2
        shapiro_val(i)=0;
    end %i
end
%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_var(Mobj, shapiro_val)
    hold on
    plot_schism_bnds(Mobj)
    colormap(jet(25))
    axis image
    box on
    auto_center
end
end

function [x3, y3, z3] = cross_product(x1,y1,z1,x2,y2,z2)

x3 = y1*z2-y2*z1;
y3 = x2*z1-x1*z2;
z3 = x1*y2-x2*y1;

end

