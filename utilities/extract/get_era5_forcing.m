function AtmForc = get_era5_forcing(Mobj, varList, src_file)
% Extract ERA5 atmopheric forcing
%
%% Syntax
% AtmForc = get_era5_forcing(Mobj, varList)
%
%% Description
% AtmForc = get_era5_forcing(Mobj, varList) extracts data from ERA5 data
%       sets based on variables in varList.
%
%% Examples
% Mobj = read_schism_hgrid(Mobj, hgrid_file);
% Mobj.force_time = (datetime(1990,12,1):hours(1):datetime(2001,2,1))';
% Mobj.force_region = Mobj.region;
% nFiles = 90;
% 
% AtmForc = get_era5_forcing(Mobj, 'prate');
% write_schism_sflux(AtmForc, 'prc', nFiles)
% 
% AtmForc = get_era5_forcing(Mobj, {'dlwrf', 'dswrf'});
% write_schism_sflux(AtmForc, 'rad', nFiles)
% 
% AtmForc = get_era5_forcing(Mobj, {'spfh', 'uwind', 'vwind', 'prmsl', 'stmp'});
% write_schism_sflux(AtmForc, 'air', nFiles)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       A datastruct containing mesh info.
% varList - variable list; cell
%       varList provides the variable names to be extracted.
% src_file - source file; char
%       the absolute filepath of ERA5 file, with **** indicating the variable name.
%       e.g., src_file = 'E:\data\ERA5_hourly_****_1980_2024.nc';
%
%% Output Arguments
% AtmForc - atmospheric forcing; datastruct
%       this datastuct contains the extracted atmospheric forcing data
%
%% Notes
% Three things to note before using this function.
%   1) All ERA5 nc files are saved by variable classification in the same directory ('ERA5_hourly_', long_name, '_1980_2024.nc', as below).
%   2) Three dimensional variables are included (lon, lat, time), and the time
%       vector is obtained from 'datenum', indicating the days since datenum(0,0,0).  
%       Note that the time vector in raw ERA5 data represents the secs
%       since 1970-01-01, and its name is 'valid_time'.
%   3) Download website (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)
% 
% ERA5_hourly_10m_u_component_of_wind_1980_2024.nc
% ERA5_hourly_10m_v_component_of_wind_1980_2024.nc   
% ERA5_hourly_2m_dewpoint_temperature_1980_2024.nc	              
% ERA5_hourly_2m_temperature_1980_2024.nc	                                   
% ERA5_hourly_mean_sea_level_pressure_1980_2024.nc	                       
% ERA5_hourly_mean_total_precipitation_rate_1980_2024.nc	                
% ERA5_hourly_surface_pressure_1980_2024.nc	                                   
% ERA5_hourly_surface_solar_radiation_downwards_1980_2024.nc	       
% ERA5_hourly_surface_thermal_radiation_downwards_1980_2024.nc	    
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024.
% Last Updated on 8 Oct 2024.
% Email: wwu@vims.edu
%
% See also: write_schism_sflux

%% Parse inputs
if ischar(varList); varList = {varList}; end
AtmForc.aimpath = Mobj.aimpath;
AtmForc.region = Mobj.force_region;
AtmForc.time = Mobj.force_time;

if nargin<3
    src_file = 'E:\ECMWF\ERA5\Hourly_single_level_chesbay\L1_data\ERA5_hourly_****_1980_2024.nc';  % download ERA5 data first!
end
%% Extract
lon_test = []; lat_test = [];

for iVar = 1:length(varList)
    varName = varList{iVar};  % input variable names

    switch lower(varName)
        case {'uwind', 'u10'} % m/s
            sflux_name = 'uwind'; % variable names in sflux nc
            long_name = '10m_u_component_of_wind';  % long variable names in ERA5 nc files (default)
            short_name = 'u10';    % short variable names in the ERA5 nc files (default)
            filepath = strrep(src_file, '****', long_name); % m/s
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);

        case {'vwind', 'v10'} % m/s
            sflux_name = 'vwind'; 
            long_name = '10m_v_component_of_wind';
            short_name = 'v10';
            filepath = strrep(src_file, '****', long_name); % m/s
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name); 

        case {'dswrf', 'swrad_down', 'short_wave_down'}  % W/m^2
            sflux_name = 'dswrf'; 
            long_name = 'surface_solar_radiation_downwards'; 
            short_name = 'ssrd'; 
            filepath = strrep(src_file, '****', long_name);  % J/m^2
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);
            varReg = varReg/3600; % convert to W/m^2

        case {'dlwrf','lwrad_down', 'long_wave_down'} % W/m^2
            sflux_name = 'dlwrf'; 
            long_name = 'surface_thermal_radiation_downwards';
            short_name = 'strd';
            filepath = strrep(src_file, '****', long_name); % J/m^2
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);
            varReg = varReg/3600; % convert to W/m^2

        case {'prmsl', 'sea_level_press', 'slp'}  % Pa, or mbar/100
            sflux_name = 'prmsl'; 
            long_name = 'mean_sea_level_pressure';
            short_name = 'msl';
            filepath = strrep(src_file, '****', long_name); % Pa, or mbar/100
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);

        case {'spfh', 'shum', 'specific_humidity'}   % specific humidity calculated at the sea level
            long_name = 'surface_pressure'; % Pa
            short_name = 'sp';
            filepath = strrep(src_file, '****', long_name);
            [~, ~, sp] = subset_era5(AtmForc, filepath, short_name);

            long_name = '2m_dewpoint_temperature';  % K
            short_name = 'd2m';
            filepath = strrep(src_file, '****', long_name);
            [lonReg, latReg, d2m] = subset_era5(AtmForc, filepath, short_name);
            sp = sp/100; d2m = d2m-273.15;
            sflux_name = 'spfh';
            varReg = calc_shum(d2m, sp);

        case {'prate', 'rain', 'precip'} 
            sflux_name = 'prate'; % kg/m^2/s
            long_name = 'mean_total_precipitation_rate'; % kg/m^2/s
            short_name = 'mtpr';
            filepath = strrep(src_file, '****', long_name);
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);

        case {'stmp',  'air_temp'} % K
            sflux_name= 'stmp';  
            long_name = '2m_temperature';
            short_name = 't2m';
            filepath = strrep(src_file, '****', long_name);
            [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name);

    end

    if sum(isnan(varReg(:)))~=0
        error('NaN values exist in the dataset, please check!')
    end

    disp(['the variable ''', sflux_name, ''' has been successfully extracted from ERA5 data set !'])
    AtmForc.(sflux_name) = varReg;
end

[AtmForc.lat, AtmForc.lon] = meshgrid(latReg, lonReg);
end

function [lonReg, latReg, varReg] = subset_era5(AtmForc, filepath, short_name)
% subset the era5 data based on your region and time period.
base_date = datenum(0,0,0); %#ok<*DATNM>

if nargin<3
    short_name = get_short_name(filepath);
end

% read dimension info
lonAll = ncread(filepath, 'lon');
latAll = ncread(filepath, 'lat');
timeAll = datetime(datevec(ncread(filepath, 'time')+base_date));

% check the time avaiability of dataset
if min(AtmForc.time) < min(timeAll) || max(AtmForc.time) > max(timeAll)
    error(['the force time exceeds the available time period in dataset (',short_name,')!'])
end

if min(AtmForc.region(1:2)) < min(lonAll) || max(AtmForc.region(1:2)) > max(lonAll) ...
        || min(AtmForc.region(3:4)) < min(latAll) || max(AtmForc.region(3:4)) > max(latAll)
    error(['the mode domain exceeds the maximum region in dataset (',short_name,')!'])
end

% check the consistency of coordinate systems
[lonNew, lon_flag] = check_lons(AtmForc.region(1:2), lonAll);
AtmForc.region(1:2) = lonNew;
if lon_flag~=0
    disp('the coordinate system of longitude is inconsisteny to your dataset!')
end

% subset era5 data
indLons = [minfind(lonAll, AtmForc.region(1)) minfind(lonAll, AtmForc.region(2))];
indLats = [minfind(latAll, AtmForc.region(3)) minfind(latAll, AtmForc.region(4))];
indTimes = [minfind(timeAll, AtmForc.time(1)) minfind(timeAll, AtmForc.time(end))];

lonReg = lonAll(min(indLons):max(indLons));
latReg = latAll(min(indLats):max(indLats));
varReg = ncread(filepath, short_name, ...
    [min(indLons) min(indLats) min(indTimes)], [abs(diff(indLons))+1 abs(diff(indLats))+1 abs(diff(indTimes))+1]); %#ok<*NASGU>

% ensure the lon/lat vectors are in ascending orders.
if ~issorted(lonReg)
    lonReg = sort(lonReg);
    varReg = flip(varReg, 1);
end
if ~issorted(latReg)
    latReg = sort(latReg);
    varReg = flip(varReg, 2);
end

% restore the longitude
switch lon_flag
    case -1
        lonReg(lonReg>180) = lonReg(lonReg>180)-360;
    case 1
        lonReg(lonReg<0) = lonReg(lonReg<0)+360;
end

end

function short_name = get_short_name(filepath)
% get the short name automatically

nc_info = ncinfo(filepath);
nc_vars = {nc_info.Variables.Name};
nc_dims = cellfun(@(x) numel(x), {nc_info.Variables.Size});
ind_var = nc_dims==max(nc_dims);
short_name = nc_vars{ind_var};

end
























