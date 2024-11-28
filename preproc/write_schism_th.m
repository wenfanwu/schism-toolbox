function write_schism_th(Mobj, prefix_name, varRaw, timeRaw)
% Write *.th file for SCHISM
%
%% Syntax 
% write_schism_th(Mobj, prefix_name, varRaw, timeRaw)
%
%% Description
% write_schism_th(Mobj, prefix_name, varRaw, timeRaw) writes a
%       <prefix_name>.th file for SCHISM
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       a datastruct containing the mesh info.
% prefix_name - file prefix name
%       prefix name of the th file. e.g., prefix_name = 'elev'.
% varRaw - raw data; double
%       2D data matrix (nt*nps) to be written.
% timeRaw - time vector; datetime
%       the corresponding time vector (datetime format).
% 
%% Notes
% this function can be used to generate elev.th, flux.th, TEM_1.th,
% SAL_1.th etc; time interval of 'timeRaw' must be greater than 'Mobj.dt' 
%
%% Author Info
% Created by Wenfan Wu,Virginia Institute of Marine Science in 2022. 
% Last Updated on 27 Nov 2024. 
% Email: wwu@vims.edu
% 
% See also: write_schism_th_nc.m

%% Parse inputs
if isnan(sum(varRaw(:)))
    error('NaNs were found in the data, please check!')
end
timeRaw = timeRaw(:);
if size(varRaw,1) ~= numel(timeRaw)
    varRaw = varRaw';
end

if size(varRaw,1) ~= numel(timeRaw)
    error('time vector is inconsistent with the data matrix!')
end

nt = size(varRaw,1);  % # of time steps
nps = size(varRaw, 2);  % # of points

time_steps = seconds(timeRaw-timeRaw(1));
if min(diff(time_steps))<Mobj.dt
    error(['the time interval should be greater than ', num2str(Mobj.dt)])
end
%% BEGIN TO WRITE
filepath = [Mobj.aimpath, prefix_name, '.th'];
fid = fopen(filepath,'wt');
format_str = repmat(['%d.', repmat('% 11.3f', [1, nps]), '\n'], [1, nt]);
varData = cat(2, time_steps(:), varRaw)';
fprintf(fid, format_str, varData(:));
fclose(fid);

disp([prefix_name, '.th has been created successfully!'])
end

