function write_schism_ptrack(Mobj, xyz_data, drop_time_list, life_day, ptrack_vars)
% Write the particle.bp file for SCHISM (Not Completed Now)
% 
%% Syntax
% write_schism_ptrack(Mobj, coords, start_time, life_day, ptrack_vars)
% 
%% Description 
% write_schism_ptrack(Mobj, coords, start_time, life_day, ptrack_vars)
% 
%% Example
% None
% 
%% Input Arguments
% ptrack_file --- the absolute filepath of 'ptrack.pth' outputed from
% 'ptrack3.exe'.
% disp_flag --- display flags. 0: no; 1: yes; Default: disp_flag = 1;
%
%% Output Arguments
% None
%
%% Notes
% The option 'life_day' is invalid for backward tracking, this is related
% to the source code of ptrack.exe in SCHISM.
% 
% This function is incomplete now, and you need to modify some code inside
% it to make sure it can work properly.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2021-11-09
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: def_schism_ptrack

%% Parse inputs
% Make sure the following parameters are the same as in 'param.nml'
dt = Mobj.dt;
nspool = 24;
ihfskip = 576;

% Check the following parameters carefully before writing
slam0 = -124; 
sfea0 = 37;
h0 = 0.01;         % min. depth for wet/dry
ndeltp = 10;      % the # of sub-divisions within each step

% Oil spill parameters
ihdf = 1; 
hdc = 3; 
horcon = 0.2;
ibuoy = 1; 
iwind = 0;
pbeach = 20;
%% Eliminate the coast nodes
cst_nodes = [Mobj.land_nodes_tot(:); Mobj.island_nodes_tot(:)];
ind_lon = Mobj.lon(cst_nodes)==xyz_data(:,1)';
ind_lat = Mobj.lat(cst_nodes)==xyz_data(:,2)';

ind_cst =  logical(sum(ind_lon & ind_lat, 1));
xyz_data(ind_cst, :) = [];

if numel(find(xyz_data(:,3)>0))~=0
    error('the release depth should be negative values!')
end
%% Parse inputs
if nargin < 4
    life_day = 10; 
    ptrack_vars = [1 0 1 0];
end
if nargin < 5
    ptrack_vars = [1 0 1 0];
end

nscreen = ptrack_vars(1);
mod_part = ptrack_vars(2);
ibf = ptrack_vars(3);
istiff = ptrack_vars(4);

switch ibf
    case 1 % forward tracking
        drop_sec_list = seconds(drop_time_list-Mobj.time(1));
        rnday = max(days(drop_time_list-Mobj.time(1)))+life_day+1; % days counting from the initiation.
    case -1 % backward tracking
        drop_sec_list = seconds(drop_time_list-Mobj.time(1));
        rnday = max(days(drop_time_list-Mobj.time(1)))+life_day+1;
end

% Coordinate option: 1: cartesian; 2: spherical
switch lower(Mobj.coord)
    case 'geographic'
        ics = 2;
    case 'cartesian'
        ics = 1;
end

%% Save the particle locations first
save([Mobj.aimpath, 'particle_locations.mat'], 'xyz_data')
disp('particle locations have been saved as MAT file!')

%% Begin to write
nps = size(xyz_data, 1);
filename = fullfile(Mobj.aimpath, 'particle.bp');

fid = fopen(filename,'w');
fprintf(fid, 'Input for ptrack*\n');
fprintf(fid, [num2str(nscreen),' !nscreen\n']);
fprintf(fid, [num2str(mod_part),' !mod_part (0: passive; 1: oil spill)\n']);
fprintf(fid, [num2str(ibf),' ! ibf (forward=1 or backward=-1 tracking)\n']);
fprintf(fid, [num2str(istiff),' ! istiff (1: fixed distance from free surface)\n']);
fprintf(fid, [num2str(ics), ' ',num2str(slam0, '%.2f'), ' ',num2str(sfea0, '%.2f'),' ! ics slam0 sfea0 (see param.nml)\n']);
fprintf(fid, [num2str(h0),' ',num2str(rnday),' ',num2str(dt),' ',num2str(nspool),' ',num2str(ihfskip),' ',num2str(ndeltp), ...
    ' ! h0,rnday,dtm,nspool,ihfskip,ndeltp (same as param.nml except for the last one, which is # of sub-divisions within each step)\n']);
fprintf(fid, [num2str(nps),'  ! # of particles\n']);
for ip = 1:nps
    switch ip
        case 1
            fprintf(fid, [num2str(ip),'   ',num2str(drop_sec_list(ip)), '  ',num2str(xyz_data(ip,1), '%.4f'),'   ',num2str(xyz_data(ip,2), '%.4f'),'   ',num2str(xyz_data(ip,3), '%.2f'), ...
                '  ! particle id, start time(s), starting x, y, and z relative to the instant surface(<=0)\n']);
        otherwise  
            fprintf(fid, [num2str(ip),'   ',num2str(drop_sec_list(ip)), '  ',num2str(xyz_data(ip,1), '%.4f'),'   ',num2str(xyz_data(ip,2), '%.4f'),'   ',num2str(xyz_data(ip,3), '%.2f'),'\n']);
    end
end
fprintf(fid, '!Oil spill parameters\n');
fprintf(fid, [num2str(ihdf),'  ',num2str(hdc, '%.1f'),'  ',num2str(horcon, '%.1f'),' ! ihdf,hdc,horcon\n']);
fprintf(fid, [num2str(ibuoy),'  ',num2str(iwind),' ! ibuoy,iwind\n']);
fprintf(fid, [num2str(pbeach, '%.1f'),' ! pbeach\n']);
fclose(fid);

disp('particle.bp has been created successfully!')

end








