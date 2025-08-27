function TideForc = add_nodal_factors(Mobj, TideForc)
% Add nodal factors for the tidal simulation
% 
%% Notes
% This function is adapted from the source code in SCHISM
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023.
% Last Updated on 24 Nov 2023.
% Email: wwu@vims.edu

% See also: ut_FUV

%% Parse inputs
load tide_fac_constants %#ok<LOAD>  % from OceanMesh2D

tideList = TideForc.tide_list;
t_s = Mobj.time(1);
t_e = Mobj.time(end);
ts = datenum(t_s); %#ok<*DATNM>
te = datenum(t_e);
tref = 0.5*(ts + te);

timeFake = linspace(ts,te,1000)';   % make a fake t with 1000 intervals
latMean = mean(Mobj.lat);
indTide = cellfun(@(x) find(contains(string(const.name), sprintf('%-4s', x))), tideList);

[F,U,V] = ut_FUV(timeFake,tref,indTide,latMean,zeros(1,4));

F = mean(F)'; % This is the average nodal factor
phs = (mean(U)+V(1,:))'*360;  % This is the average nodal correction + astronomical argument at beginning of simulation

% Make sure phase between 0 and 360 for aesthetic purposes
while any(phs < 0)
    phs(phs<0) = phs(phs<0) + 360;
end
TideForc.nodal_factor = F;
TideForc.eq_arg = phs;

end
