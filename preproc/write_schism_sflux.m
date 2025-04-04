function write_schism_sflux(AtmForc, suffix_name, time_steps, start_idx, ind_dst)
% Write sflux nc files for SCHISM
%
%% Syntax 
% write_schism_sflux(AtmForc, suffix_name)
% write_schism_sflux(AtmForc, suffix_name, time_steps)
% write_schism_sflux(AtmForc, suffix_name, time_steps, start_idx)
% write_schism_sflux(AtmForc, suffix_name, time_steps, start_idx, ind_dst)
%
%% Description
% write_schism_sflux(AtmForc, suffix_name) writes the sflux files in NetCDF format
% write_schism_sflux(AtmForc, suffix_name, time_steps) specifies the time steps in each file
% write_schism_sflux(AtmForc, suffix_name, time_steps, start_idx) specifies the start index of the files
% write_schism_sflux(AtmForc, suffix_name, time_steps, start_idx, ind_dst) specifies the index of dataset.
% 
%% Example
% write_schism_sflux(AtmForc, 'rad')
% write_schism_sflux(AtmForc, 'rad', 900)
% write_schism_sflux(AtmForc, 'rad', 900, 1)
% 
%% Input Arguments
% AtmForc - atmospheric forcing; datastruct
%       this datastruct stores atmospheric forcing data.
% suffix_name - filename suffix; char
%       the suffix filename of NetCDF files (air/prc/rad).
% time_steps - the # of time steps (optional); numeric
%       the # of time steps in each NetCDF file. time_steps should be less
%       than 1000, which is the requirement in SCHISM code by default.
%       Default: time_steps = 995. The upper limit of time_steps is set to
%       995 rather than 1000, in order to make sure the time steps of the
%       last file will not exceed 1000.
% start_idx - the starrt index (optional); numeric
%       the start index of the first stacked NetCDF files. default: start_idx = 1
% ind_dst - the index of dataset (optional); numeric
%       the index of used dataset. Default: ind_dst = 1.  
%       In SCHISM, if multiple datasets are available in sflux files, the
%       files with larger index with overwrite the previous one when they
%       are overlapped in space. 
% 
%% Output Arguments
% None
%
%% Notes
% If the simulation period is long and the model domain is extremely large,
% you may encounter memory issues when attempting to extract the entire
% dataset at once. In this case, it is recommended to write the sflux files
% in multiple continuous time windows by incrementing the "start_idx".
% 
% Change 'max_files' and 'max_times' in 'sfux_9c.F90' if needed.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 3 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: get_era5_forcing

%% Parse inputs
if nargin < 3; time_steps = 995; end
if nargin < 4; start_idx = 1; end
if nargin < 5; ind_dst = 1; end
if time_steps > 995; error('time steps should not exceed 1000!'); end
if time_steps == 1; error('time steps should be greater than 1!'); end

%% Check dimensions and NaNs
field_list = fieldnames(AtmForc);
ind_vars = ~contains(field_list, {'lon', 'lat', 'time', 'aimpath', 'region'});
varList = field_list(ind_vars);

% check dimensions
if dimnum(AtmForc.lon)==1 || dimnum(AtmForc.lat)==1 || ...
        numel(unique(AtmForc.lon(1,:)))~=1 || numel(unique(AtmForc.lat(:,1)))~=1
    error('lon/lat should be created by MESHGRID function!')
end
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

sflux_path = fullfile(AtmForc.aimpath, 'sflux\');
if exist(sflux_path,'dir')~=7
    mkdir(sflux_path)
end
%% Check time steps
[nLons, nLats] = size(AtmForc.lon);
total_steps = length(AtmForc.time);
nFiles = ceil(total_steps/time_steps);  % save the last file (without enough time steps) seperately by default.

if mod(total_steps, time_steps) == 1  % it will cause the last file have no time dimension.
    nFiles = nFiles-1;
end

% change the two parameters in 'sflux_9c.F90' if necessary
max_files = 9999;             % max. total # of nc files  (default: 9999).
max_times = 100000;        % max. # of time records from all files (default: 100000).

warning on
if nFiles > max_files; warning('nFiles is greater than the max # of nc files!'); end
if total_steps > max_times; warning('nFiles is greater than max. # of time records from all files!'); end

n_digits = length(num2str(floor(abs(nFiles))));
%% sflux nc files
dtype = 'single'; nc_fmt = 'classic';  % define the NetCDF format

for iFile = 1:nFiles
    progressbar(iFile/nFiles)
    
    beg_ind = 1+time_steps*(iFile-1);
    end_ind = time_steps*iFile;
    if iFile == nFiles; end_ind = total_steps; end

    sub_time = AtmForc.time(beg_ind:end_ind);
    base_date = dateshift(sub_time(1), 'start', 'days');  % the 'hour' component of 'base_date' is unused.
    time = datenum(sub_time)-datenum(base_date); %#ok<*DATNM>
    
    nTimes = length(time);
    base_date_str = datestr(base_date, 'yyyy-mm-dd'); %#ok<*DATST>
    
    fileName = [sflux_path, 'sflux_', suffix_name, '_', num2str(ind_dst), '.', num2str(start_idx-1+iFile, ['%0',num2str(max(n_digits,4)),'d']),'.nc'];
    if exist(fileName,'file')==2; delete(fileName); end
    
    %---------TIME PART
    nccreate(fileName,'time','Dimensions',{'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
    ncwriteatt(fileName,'time','long_name','Time');
    ncwriteatt(fileName,'time','standard_name','time');
    ncwriteatt(fileName,'time','units', ['days since ', base_date_str]);
    ncwriteatt(fileName,'time','base_date', int32([year(base_date) month(base_date) day(base_date) hour(base_date)]));
    ncwrite(fileName,'time', time);
    
    nccreate(fileName,'lon','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats},'Datatype', dtype,'Format', nc_fmt)
    ncwriteatt(fileName,'lon','long_name','Longitude');
    ncwriteatt(fileName,'lon','standard_name','longitude');
    ncwriteatt(fileName,'lon','units', 'degree_east');
    ncwrite(fileName,'lon', AtmForc.lon);
    
    nccreate(fileName,'lat','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats},'Datatype', dtype,'Format', nc_fmt)
    ncwriteatt(fileName,'lat','long_name','Latitude');
    ncwriteatt(fileName,'lat','standard_name','latitude');
    ncwriteatt(fileName,'lat','units', 'degree_north');
    ncwrite(fileName,'lat', AtmForc.lat);
    
    ncwriteatt(fileName,'/','Conventions', 'CF-1.0');
    
    %---------VARIABLE PART
    % ------- AIR
    if strcmp(suffix_name, 'air')
        nccreate(fileName,'uwind','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'uwind','long_name','Surface Eastward Air Velocity (10m AGL)');
        ncwriteatt(fileName,'uwind','standard_name','eastward_wind');
        ncwriteatt(fileName,'uwind','units','m/s');
        ncwrite(fileName,'uwind', AtmForc.uwind(:,:,beg_ind:end_ind));
        
        nccreate(fileName,'vwind','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'vwind','long_name','Surface Northward Air Velocity (10m AGL)');
        ncwriteatt(fileName,'vwind','standard_name','northward_wind');
        ncwriteatt(fileName,'uwind','units','m/s');
        ncwrite(fileName,'vwind', AtmForc.vwind(:,:,beg_ind:end_ind));
        
        nccreate(fileName,'prmsl','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'prmsl','long_name','Pressure reduced to MSL');
        ncwriteatt(fileName,'prmsl','standard_name','air_pressure_at_sea_level');
        ncwriteatt(fileName,'prmsl','units','Pa');
        ncwrite(fileName,'prmsl', AtmForc.prmsl(:,:,beg_ind:end_ind));
        
        nccreate(fileName,'stmp','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'stmp','long_name','Surface Air Temperature (2m AGL)');
        ncwriteatt(fileName,'stmp','standard_name','air_temperature');
        ncwriteatt(fileName,'stmp','units','K');
        ncwrite(fileName,'stmp', AtmForc.stmp(:,:,beg_ind:end_ind));
        
        nccreate(fileName,'spfh','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'spfh','long_name','Surface Specific Humidity (2m AGL)');
        ncwriteatt(fileName,'spfh','standard_name','specific_humidity');
        ncwriteatt(fileName,'spfh','units','1');
        ncwrite(fileName,'spfh', AtmForc.spfh(:,:,beg_ind:end_ind));
    end
    % ------- PRC
    if strcmp(suffix_name, 'prc')
        nccreate(fileName,'prate','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'prate','long_name','Surface Precipitation Rate');
        ncwriteatt(fileName,'prate','standard_name','precipitation_flux');
        ncwriteatt(fileName,'prate','units','kg/m^2/s');
        ncwrite(fileName,'prate', AtmForc.prate(:,:,beg_ind:end_ind));
    end
    % ------- RAD
    if strcmp(suffix_name, 'rad')
        nccreate(fileName,'dlwrf','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'dlwrf','long_name','Downward Long Wave Radiation Flux');
        ncwriteatt(fileName,'dlwrf','standard_name','surface_downwelling_longwave_flux_in_air');
        ncwriteatt(fileName,'dlwrf','units','W/m^2');
        ncwrite(fileName,'dlwrf', AtmForc.dlwrf(:,:,beg_ind:end_ind));
        
        nccreate(fileName,'dswrf','Dimensions',{'nx_grid', nLons, 'ny_grid', nLats, 'time', nTimes},'Datatype', dtype,'Format', nc_fmt)
        ncwriteatt(fileName,'dswrf','long_name','Downward Short Wave Radiation Flux');
        ncwriteatt(fileName,'dswrf','standard_name','surface_downwelling_shortwave_flux_in_air');
        ncwriteatt(fileName,'dswrf','units','W/m^2');
        ncwrite(fileName,'dswrf', AtmForc.dswrf(:,:,beg_ind:end_ind));
    end
end

disp([num2str(nFiles), ' stacked sflux files has been created !'])
%% sflux_inputs.txt (default)
fid = fopen([sflux_path, 'sflux_inputs.txt'],'wt');
fprintf(fid, '&sflux_inputs\n');   
fprintf(fid, '/');   
fclose(fid);

if isfield(AtmForc, 'dataset')
    dst = AtmForc.dataset;
else
    dst = 'unknown';
end
fid = fopen([sflux_path, 'README.md'],'wt');
fprintf(fid, ['product name: ', dst, '\n']);  % TBD
fprintf(fid, ['timespan: ',datestr(AtmForc.time(1), 'yyyymmddTHHMMSSZ'),' to ', ...
    datestr(AtmForc.time(end), 'yyyymmddTHHMMSSZ'), '\n']);  
fprintf(fid, ['region: ',num2str(AtmForc.region(1)), '-',num2str(AtmForc.region(2)), '; ', ...
    num2str(AtmForc.region(3)), '-', num2str(AtmForc.region(4)), '\n']);  
fclose(fid);
end