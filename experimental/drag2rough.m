function rough_p = drag2rough(Mobj, Cd, dzb_min)
% Convert from Cd to roughness [m] in SCHISM model

if nargin < 3
    dzb_min = 0.5; % min. bottom boundary layer thickness [m] in param.nml
end

bthick = get_schout_btm(Mobj, Mobj.depLayers,1)+Mobj.depth';  % bottom layer thickness [m]
bthick = max(dzb_min, bthick);
rough_p = bthick(:).*exp(-0.4./sqrt(Cd));

end