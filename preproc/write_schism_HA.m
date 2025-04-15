function write_schism_HA(Mobj, HA_tides, HA_dvars, HA_flags)
% Write the harm.in file for SCHISM
%
%% Syntax
% write_schism_HA(Mobj, HA_tides, HA_dvars) 
% write_schism_HA(Mobj, HA_tides, HA_dvars, HA_flags)
%
%% Description
% write_schism_HA(Mobj, HA_tides, HA_dvars) writes the harm.in file for SCHISM
% write_schism_HA(Mobj, HA_tides, HA_dvars, HA_flags) specifies the flags for elev and velocity
%
%% Examples 
% HA_tides = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
% HA_dvars = [10 40 24 0];
% HA_flags = [1 0];
% write_schism_HA(Mobj, HA_tides, HA_dvars, HA_flags)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% HA_tides - tide constituents; cell
%       the tide constituents used for harmonic analysis;
% HA_dvars - date variables; numeric
%       the date variables used for harmonic analysis; 
%       HA_dvars = [start_day, end_day, stride_step, de_tend];
% HA_flags - ha flags; numeric
%       the flags used for elevation and velocity. 
%       Default: HA_flags = [1, 0]; (latter not active).
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 15 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_bctides

%% Parse inputs
if nargin < 3; HA_dvars = [10 40 24 0]; end
if nargin < 4; HA_flags = [1 0]; end

start_day = HA_dvars(1); end_day = HA_dvars(2);
stride_step = HA_dvars(3); de_tend = HA_dvars(4);
elev_flag = HA_flags(1); vel_flag = HA_flags(2);

%% Tidal parameters
load('tide_fac_constants.mat', 'const');
tide_pool = cellstr(strtrim(string(const.name)));

nTides = length(HA_tides);
tide_freqs = nan(nTides,1);
for iTide = 1:nTides
    tide_name = strtrim(HA_tides{iTide});
    ind_tide = find(strcmpi(tide_pool, tide_name));

    tide_freqs(iTide) = const.freq(ind_tide)/3600*2*pi;  % rad/s
end

TideForc.tide_list = HA_tides;
TideForc = add_nodal_factors(Mobj, TideForc);

nodal_factor = TideForc.nodal_factor;
eq_arg = TideForc.eq_arg;

%% Begin to write
filepath = [Mobj.aimpath, 'harm.in'];
fid = fopen(filepath,'w');

fprintf(fid, '%d !# of freq\n', nTides);
for iTide = 1:nTides
    % tidal constituent name
    fprintf(fid, [upper(HA_tides{iTide}), '\n']);

    % angular frequency, nodal factor, earth equilibrium argument (in degrees);
    if iTide == 1
        fprintf(fid, '%12.6e% 8.5f% 9.5f !angular freq, nodal factor, argument (deg)\n', tide_freqs(iTide), nodal_factor(iTide), eq_arg(iTide));
    else
        fprintf(fid, '%12.6e% 8.5f% 9.5f\n', tide_freqs(iTide), nodal_factor(iTide), eq_arg(iTide));
    end
end
fprintf(fid, '%d %d %d %.2f !starting day, end day, stride in step, de-tending\n', start_day, end_day, stride_step, de_tend);
fprintf(fid, '%d %d !do HA for elev, vel (latter not active)\n', elev_flag, vel_flag);

fclose(fid);
disp('harm.in has been created successfully!')

end