function flux_flags = def_schism_fluxflag(Mobj, nRegs)
% Define fluxflags on the SCHISM grid
%
%% Syntax
% flux_flags = def_schism_fluxflag(Mobj)
% flux_flags = def_schism_fluxflag(Mobj, nRegs)
%
%% Description
% flux_flags = def_schism_fluxflag(Mobj) defines fluxflags
% flux_flags = def_schism_fluxflag(Mobj, nRegs) specifies the # of defined regions
%
%% Examples
% figure
% disp_schism_hgrid(Mobj, [0 0])
% flux_flags = def_schism_fluxflag(Mobj, 2);
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% nRegs - the # of regions; double
%       the # of regions to be defined. default: nReg = 1;
%
%% Output Arguments
% flux_flags - flux flags; double
%       the flux flags (nElems*1) used to calculate transect flux.
%
%% Author Info.
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 18 Dec 2024.
% Email: wwu@vims.edu
% 
% See also: def_schism_mask

%% Parse inputs
if nargin < 2
    nRegs = 1;
end
num_list = 1:nRegs; % can be changed
nRegs = numel(num_list);
%% Define
flux_flags = (min(num_list)-1)*ones(Mobj.nElems, 1);

hold on
for iReg = 1:nRegs
    disp('draw a polygon on the map and press ENTER')
    geo_handle = drawpolygon;
    lonRoi = geo_handle.Position(:,1)';
    latRoi = geo_handle.Position(:,2)';
    delete(geo_handle)
    
    ind_flag = inpolygon(Mobj.lonc, Mobj.latc, lonRoi, latRoi);
    lon_list = Mobj.lonc(ind_flag);
    lat_list = Mobj.latc(ind_flag);
    scatter(lon_list, lat_list, 5, 'magenta', 'filled')
    flux_flags(ind_flag) = num_list(iReg);
end
end