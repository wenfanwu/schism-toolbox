function write_schism_sflux(AtmForc, suffix_name, nFiles, begNum, nDataset)
% Write sflux nc files for SCHISM
%
%% Syntax 
% write_schism_sflux(AtmForc, suffix_name, nFiles)
% write_schism_sflux(AtmForc, suffix_name, nFiles, begNum)
% write_schism_sflux(AtmForc, suffix_name, nFiles, begNum, nDataset)
%
%% Description
% write_schism_sflux(AtmForc, suffix_name, nFiles) writes the sflux
%       NetCDF files for SCHISM model.
% write_schism_sflux(AtmForc, suffix_name, nFiles, begNum) specifies the
%       beginning index for the NetCDF files.
% write_schism_sflux(AtmForc, suffix_name, nFiles, begNum, nDataset)
%       specifies the numbering of dataset.
% 
%% Example
% write_schism_sflux(AtmForc, 'rad', 30)
% 
%% Input Arguments
% AtmForc - atmospheric forcing; datastruct
%       this datastruct stores atmospheric forcing data.
% suffix_name - filename suffix; char
%       used to specify the suffix name of NetCDF files (air/prc/rad).
% nFiles - # of files; double
%       the # of stacked files. Since the 'hour' component of
%       base_date is unused in each stacked file, thus this function will
%       slightly adjuct the nFiles to ensures that each file starts at 0 o 'clock
%       in the day.
% 
%% Output Arguments
% None
%
%% Notes
% If the simulation period is too long, it is recommended to prepare the
% netcdf files in multiple temporal segments, just speficy "depNum"
% greater than 1.
% 
% Change 'max_files' and 'max_times' in 'sfux_9c.F90' if needed.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 24 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: get_era5_forcing

%% Parse inputs
if nargin < 4
    begNum = 1;
end
if nargin < 5
    nDataset = 1;
end
%% Check dimensions and NaNs
fieldList = fieldnames(AtmForc);
ind_vars = ~contains(fieldList, {'lon', 'lat', 'time', 'aimpath', 'region'});
varList = fieldList(ind_vars);

% check dimensions
if numel(find(isnan(AtmForc.lon(:))))~=0
    error('NaN values exist in the longitude matrix!')
end
if numel(find(isnan(AtmForc.lat(:))))~=0
    error('NaN values exist in the latitude matrix!')
end
if numel(find(isnat(AtmForc.time(:))))~=0
    error('NaN values exist in the time vector!')
end

% check variables
nVars = length(varList);
for iVar = 1:nVars
    varName = varList{iVar};
    varData = AtmForc.(varName);
    if sum(size(varData)==[size(AtmForc.lat) length(AtmForc.time)])~=3
        error(['lon/lat/time dimension is inconsistent with the variable (', varName, ') !'])
    end
    if numel(find(isnan(varData(:))))~=0
        error(['NaN values exist in the variable (', varName, ') !'])
    end
end

% ensure the time interval is even
time_interval = unique(diff(AtmForc.time));
if numel(time_interval) ~= 1
    error('the time interval for the atmospheric forcing is uneven!')
end
time_cell = days(1)/time_interval;

sflux_path = fullfile(AtmForc.aimpath, 'sflux\');
if exist(sflux_path,'dir')~=7
    mkdir(sflux_path)
end
%% Adjust nFiles
nLons = size(AtmForc.lon,1);
nLats = size(AtmForc.lon,2);
total_times = length(AtmForc.time);

time_steps = floor(1000/time_cell)*time_cell;
nFiles_min = ceil(total_times/time_steps);

if nargin < 3
    nFiles = nFiles_min;
    disp(['nFiles is set to ', num2str(nFiles), ' automatically'])
else
    time_steps = floor(total_times/nFiles);
    if time_steps<time_cell 
        warning on
        warning('nFiles is too large to ensure each file starts at zero hour')
        time_steps = time_cell;
    end
    time_steps = floor(time_steps/time_cell)*time_cell;

    nFiles = ceil(total_times/time_steps);
    disp(['nFiles is adjusted to ', num2str(nFiles)])
end

if time_steps > 1000
    error(['nFiles should be less than ', num2str(nFiles_min),' since time records in each nc files cannot exceed 1000 by default'])
end

digits_num = length(num2str(floor(abs(nFiles))));
%% Check data entries
% Change theses two parameters in 'sfux_9c.F90' if needed.
max_files = 9999;             % max. total # of nc files  (default: 9999).
max_times = 300000;        % max. # of time records from all files (default: 100000).

if nFiles>max_files
    warning('nFiles is greater than the max # of nc files!')
end
if total_times>max_times
    warning('nFiles is greater than max. # of time records from all files!')
end
%% sflux nc files
for iFile = 1:nFiles
    begind = 1+time_steps*(iFile-1);
    endind = time_steps*iFile;
    if iFile == nFiles
        endind = total_times;
    end
    sub_time = AtmForc.time(begind:endind);
    base_date = dateshift(sub_time(1), 'start', 'days');
    time = (datenum(sub_time)-datenum(base_date))'; %#ok<*DATNM>
    
    nTimes = length(time);
    begtime_str = datestr(sub_time(1), 'yyyy-mm-dd'); %#ok<*DATST>
    
    fileName = [sflux_path, 'sflux_', suffix_name, '_',num2str(nDataset),'.',num2str(begNum-1+iFile, ['%0',num2str(max(digits_num,4)),'d']),'.nc'];
    if exist(fileName,'file')==2; delete(fileName); end
    
    %---------TIME PART
    nccreate(fileName,'time','Dimensions',{'time', nTimes},'Datatype','single','Format','classic')
    ncwriteatt(fileName,'time','long_name','Time');
    ncwriteatt(fileName,'time','standard_name','time');
    ncwriteatt(fileName,'time','units', ['days since ', begtime_str]);
    ncwriteatt(fileName,'time','base_date', int32([year(base_date) month(base_date) day(base_date) hour(base_date)]));
    ncwrite(fileName,'time', time);
    
    nccreate(fileName,'lon','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats},'Datatype','single','Format','classic')
    ncwriteatt(fileName,'lon','long_name','Longitude');
    ncwriteatt(fileName,'lon','standard_name','longitude');
    ncwriteatt(fileName,'lon','units', 'degree_east');
    ncwrite(fileName,'lon', AtmForc.lon);
    
    nccreate(fileName,'lat','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats},'Datatype','single','Format','classic')
    ncwriteatt(fileName,'lat','long_name','Latitude');
    ncwriteatt(fileName,'lat','standard_name','latitude');
    ncwriteatt(fileName,'lat','units', 'degree_north');
    ncwrite(fileName,'lat', AtmForc.lat);
    
    ncwriteatt(fileName,'/','Conventions', 'CF-1.0');
    
    %---------VARIABLE PART
    % ------- AIR
    if strcmp(suffix_name, 'air')
        nccreate(fileName,'uwind','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'uwind','long_name','Surface Eastward Air Velocity (10m AGL)');
        ncwriteatt(fileName,'uwind','standard_name','eastward_wind');
        ncwriteatt(fileName,'uwind','units','m/s');
        ncwrite(fileName,'uwind', AtmForc.uwind(:,:,begind:endind));
        
        nccreate(fileName,'vwind','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'vwind','long_name','Surface Northward Air Velocity (10m AGL)');
        ncwriteatt(fileName,'vwind','standard_name','northward_wind');
        ncwriteatt(fileName,'uwind','units','m/s');
        ncwrite(fileName,'vwind', AtmForc.vwind(:,:,begind:endind));
        
        nccreate(fileName,'prmsl','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'prmsl','long_name','Pressure reduced to MSL');
        ncwriteatt(fileName,'prmsl','standard_name','air_pressure_at_sea_level');
        ncwriteatt(fileName,'prmsl','units','Pa');
        ncwrite(fileName,'prmsl', AtmForc.prmsl(:,:,begind:endind));
        
        nccreate(fileName,'stmp','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'stmp','long_name','Surface Air Temperature (2m AGL)');
        ncwriteatt(fileName,'stmp','standard_name','air_temperature');
        ncwriteatt(fileName,'stmp','units','K');
        ncwrite(fileName,'stmp', AtmForc.stmp(:,:,begind:endind));
        
        nccreate(fileName,'spfh','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'spfh','long_name','Surface Specific Humidity (2m AGL)');
        ncwriteatt(fileName,'spfh','standard_name','specific_humidity');
        ncwriteatt(fileName,'spfh','units','1');
        ncwrite(fileName,'spfh', AtmForc.spfh(:,:,begind:endind));
    end
    % ------- PRC
    if strcmp(suffix_name, 'prc')
        nccreate(fileName,'prate','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'prate','long_name','Surface Precipitation Rate');
        ncwriteatt(fileName,'prate','standard_name','precipitation_flux');
        ncwriteatt(fileName,'prate','units','kg/m^2/s');
        ncwrite(fileName,'prate', AtmForc.prate(:,:,begind:endind));
    end
    % ------- RAD
    if strcmp(suffix_name, 'rad')
        nccreate(fileName,'dlwrf','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'dlwrf','long_name','Downward Long Wave Radiation Flux');
        ncwriteatt(fileName,'dlwrf','standard_name','surface_downwelling_longwave_flux_in_air');
        ncwriteatt(fileName,'dlwrf','units','W/m^2');
        ncwrite(fileName,'dlwrf', AtmForc.dlwrf(:,:,begind:endind));
        
        nccreate(fileName,'dswrf','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype','single','Format','classic')
        ncwriteatt(fileName,'dswrf','long_name','Downward Short Wave Radiation Flux');
        ncwriteatt(fileName,'dswrf','standard_name','surface_downwelling_shortwave_flux_in_air');
        ncwriteatt(fileName,'dswrf','units','W/m^2');
        ncwrite(fileName,'dswrf', AtmForc.dswrf(:,:,begind:endind));
    end
end

disp([num2str(nFiles), ' stacked sflux files has been created !'])
%% sflux_inputs.txt (default)
fid = fopen([sflux_path, 'sflux_inputs.txt'],'wt');
fprintf(fid, '&sflux_inputs\n');   
fprintf(fid, '/');   
fclose(fid);

fid = fopen([sflux_path, 'README.md'],'wt');
fprintf(fid, 'product name: ERA5\n');  % TBD
fprintf(fid, ['timespan: ',datestr(AtmForc.time(1), 'yyyymmddTHHMMSSZ'),' to ', ...
    datestr(AtmForc.time(end), 'yyyymmddTHHMMSSZ'), '\n']);  
fprintf(fid, ['region: ',num2str(AtmForc.region(1)), '-',num2str(AtmForc.region(2)), '; ', ...
    num2str(AtmForc.region(3)), '-', num2str(AtmForc.region(4)), '\n']);  
fclose(fid);
end
