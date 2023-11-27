function flux_flags = def_schism_fluxflag(Mobj, num_regs)
%
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
% Mobj --- Mesh object
% num_regs --- the # of defined regions.
%
%% Output Arguments
% flux_flags --- flux flags@elements
%% Notes
%
%
%% Author Info.
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2022-10-02.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
num_list = 1:num_regs; % consecutive positive integers
flux_flags = (min(num_list)-1)*ones(Mobj.nElems, 1);

for N = 1:num_regs
    figure('Color', 'w')
    disp_schism_hgrid(Mobj)
    axis image
    colorbar;
    colormap(jet(25))
    
    disp('draw a polygon on the map to select the region')
    geo_handle = drawpolygon;
    lonRoi = geo_handle.Position(:,1)';
    latRoi = geo_handle.Position(:,2)';
    close gcf
    
    ind_flag = inpolygon(Mobj.lonc, Mobj.latc, lonRoi, latRoi);
    lon_list = Mobj.lonc(ind_flag);
    lat_list = Mobj.latc(ind_flag);
    ind_elems = geomin(Mobj.lonc, Mobj.latc, lon_list, lat_list);
    flux_flags(ind_elems) = num_list(N);
    
    close gcf
end

%% Display
figure('Color', 'w')
disp_schism_var(Mobj, flux_flags)
caxis([min(num_list)-0.5 max(flux_flags)+0.5])
axis image
end