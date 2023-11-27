function write_schism_station_in(Mobj, xyz_data, switch_flags)
% Write station.in file for SCHISM model
% 
%% Syntax
% write_schism_station_in(Mobj, xyz_data, switch_flags)
%
%% Description 
% write_schism_station_in(Mobj, xyz_data, switch_flags) creates the
% station.in file for SCHISM.
%
%% Examples
%  write_schism_station_in(Mobj, xyz_data, switch_flags)
%
%% Input Arguments
% Mobj --- the mesh object
% xyz_data --- a N*3 matrix contains the coordinates of station points.
% xyz_data = [xsta(:), ysta(:), zsta(:)].
% swtich_flags --- flags used to determine which variables will be output.
% Default: switch_flags = [1 1 1 1 1 1 1 1 1];
%
%% Output Arguments
% None
% 
%% Notes
%  Format: station #,x,y,z; 
%  if ics=2, x,y are degrees in lon/lat. z is z-coord (not distance from surface!). 
%  For 3D variables, code will extrapolate above surface/below bottom if necessary.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-26.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if nargin < 3
    switch_flags = [1 1 1 1 1 1 1 1 1]; % on (1) | off(0) flags for elev, air pressure, windx, windy, T, S, u, v, w
end
if numel(switch_flags) ~=9
    default_flags = zeros(1,9);
    default_flags(1:length(switch_flags)) = switch_flags;
    switch_flags = default_flags;
    warning('Specified switch is insufficient, the unspecified will not be output by default')
end

nsta = size(xyz_data,1);  % # of stations
xsta = xyz_data(:,1);
ysta = xyz_data(:,2);
zsta = xyz_data(:,3);  % (from vertical datum; <0 is below)

filepath = [Mobj.aimpath, 'station.in'];
fid = fopen(filepath,'wt');
fprintf(fid, '%d% d% d% d% d% d% d% d% d !on (1)|off(0) flags for elev, air pressure, windx, windy, T, S, u, v, w\n', switch_flags);
fprintf(fid, '%d\n', nsta);
for ii = 1:nsta
    fprintf(fid, '%d% 11.6f% 8.5f% d\n', ii, xsta(ii), ysta(ii), zsta(ii));
end
fclose(fid);

end

% NOTES
% 1 1 1 1 1 1 1 1 1 !on (1)|off(0) flags for elev, air pressure, windx, windy, T, S, u, v, w
% nsta ! # of stations
% do i=1,np
%     i,xsta(i),ysta(i),zsta(i) !zsta(i) is z-coordinates (from vertical datum; <0 is below)
% enddo    









