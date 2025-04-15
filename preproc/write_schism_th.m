function write_schism_th(Mobj, prefix_name, varRaw, timeRaw)
% Write the *.th file for SCHISM
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
%       the datastruct used to store mesh info.
% prefix_name - fprefix name; char
%       prefix name of the *.th file. e.g., prefix_name = 'elev' or
%       'elev.th'; the file extension will be omitted automatically.
% varRaw - raw data; numeric
%       the 2D data matrix (nt*nps) to be written.
% timeRaw - time vector; datetime
%       the corresponding time vector, the time interval of 'timeRaw' must
%       be greater than 'Mobj.dt'.
% 
%% Notes
% This function can be used to generate elev.th, flux.th, TEM_1.th,
% SAL_1.th etc
%
%% Author Info
% Created by Wenfan Wu,Virginia Institute of Marine Science in 2022. 
% Last Updated on 5 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_th_nc

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.th$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '.th']);

%% Check inputs
timeRaw = timeRaw(:);
if any(isnan(varRaw(:))); error('NaNs were found in the input data!'); end
if size(varRaw,1) ~= numel(timeRaw); varRaw = varRaw'; end
if size(varRaw,1) ~= numel(timeRaw); error('time vector is inconsistent with the data matrix!'); end

[~, nps] = size(varRaw);  % # of time steps
time_secs = seconds(timeRaw-timeRaw(1));
if min(diff(time_secs))<Mobj.dt; error(['time interval should be greater than ', num2str(Mobj.dt)]); end

%% BEGIN TO WRITE
fid = fopen(filepath,'w');
fprintf(fid, ['%d.', repmat('% 11.3f', [1, nps]), '\n'], cat(2, time_secs(:), varRaw)');
fclose(fid);

disp([prefix_name, '.th has been created successfully!'])
end

