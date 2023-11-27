function D = get_hycom_online(aimpath, region, timeTick, varList, URL)
% Download the HYCOM data in a flexible way
%
%% Syntax
% D = get_hycom_online(aimpath, region, timeTick)
% D = get_hycom_online(aimpath, region, timeTick, varList)
% D = get_hycom_online(aimpath, region, timeTick, varList, URL)
%
%% Description
% D = get_hycom_online(aimpath, region, timeTick) downloads HYCOM data for
% a particular moment and region into a specified folder.
%
% D = get_hycom_online(aimpath, region, timeTick, varList) specifies the
% required variables.
%
% D = get_hycom_online(aimpath, region, timeTick, varList, URL) specifies
% the HYCOM product.
%
%% Example-1: Download HYCOM data of a particular moment
% clc;clearvars
% aimpath = 'E:/data/';
% region = [190 240 -5 5]; % Nino3.4
% timeTick = datetime(2010,1,1);
% varList = {'ssh','temp','salt','uvel','vvel'};
% D = get_hycom_online(aimpath, region, timeTick, varList);
%
%% Example-2: Download HYCOM data in bulk
% clc;clearvars
% aimpath = 'E:/data/';
% region = [117.5 122.5 37 41]; % the Bohai Sea
% timeList = datetime(2020,1,1):hours(3):datetime(2020,2,1);
% varList = {'ssh','temp','salt','uvel','vvel'};
%
% nTimes = numel(timeList);
% for iTime = 1:nTimes
%     timeTick = timeList(iTime);
%     D = get_hycom_online(aimpath, region, timeTick, varList);
% end
%
%% Example-3: Download data from a specified HYCOM product
% clc;clearvars
% aimpath = 'E:/data/';
% region = [261 280 17.5 32.5]; % the Gulf of Mexico
% timeTick = datetime(2010,1,1);
% varList = {'ssh','temp','salt','u','v'};
% URL = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?';
% D = get_hycom_online(aimpath, region, timeTick, varList, URL);
%
%% Input Arguments
% aimpath --- the directory where the HYCOM data is stored. It doesn't
% matter if this directory name ends with a backslash.
%
% region --- the region of interest. e.g. region = [lon_south, lon_north, lat_south, lat_north];
% the longitude should be in [0, 360], while latitude is in [-80 80].
%
% timeTick --- the specified time with datetime format. e.g. timeTick = datetime(2010,1,1);
%
% varList --- variable list. Default: varList = {'ssh','temp','salt','u','v'};
%
%% Output Arguments
% D --- a datastruct containing all the variables you need. Note that the
% dev_time field means the deviation (in hours) between the actual time of
% the downloaded data and your specified time.
%
%% Notes
% There are three things to note before you use this function:
%
% (1) this function aims to download the HYCOM data of a particular moment,
% and it will search the HYCOM data at the nearest moment relative to your
% given time.
%
% (2) before downloading, this function will check if you have ever made
% the same request, if so, it will load the available one directly.
%
% (3) this function integrates 13 types of HYCOM products, all of which
% have latitude vectors from -80 to 80. However, there are 8 products
% whose longitude vectors are from 0 to 360, whereas those of the other 5
% products are from -180 to 180 (see the bottom of this function). This
% inconsistency slightly hinders our data reading, especially when your
% given spatial region intersects the prime meridian.
%
%% Tips
% The network of HYCOM website seems to be unstable, so it may take
% a long time to download data sometimes, or encountered netCDF errors,
% just re-run the function in this case.
%
% A possible method to accelerate this function is to save the dimension
% info (lon, lat, time, depth) of different HYCOM products as MAT files in
% advance, but it may reduce the conciseness of this function.
%
% What if the HYCOM data of the particular time is missing? You can increase
% the 'tole_time' parameter at Line 175 to replace missing data with data on
% adjacent days. It means tolerance time bias. e.g. tole_time = days(3);
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021.
% Last Updated on 2022-05-20.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also: ncread

%% Parse inputs
tic
if exist(aimpath,'dir')~=7
    disp('the aimpath does not exist and has been created automatically')
    mkdir(aimpath);
end

if min(region(1:2))<0 || max(region(1:2))>360 || min(region(3:4))<-80 || max(region(3:4))>80
    error('The longitudes should be in [0, 360], while latitudes is in [-80 80].')
end

varBase = {'ssh','temp','salt','u','v'};
stdBase = {'surf_el', 'water_temp','salinity','water_u','water_v'};

if nargin<4
    varList = {'ssh','temp','salt','u','v'};                                              % variable name list
    stdList = {'surf_el', 'water_temp','salinity','water_u','water_v'};  % standard name list
end

varList = lower(varList);
ind_vars = cellfun(@(x) find(contains(varBase, x)), varList);
if numel(ind_vars)~=numel(varList)
    warning on
    warning('Some variable names are unrecognized!')
end
varList = varBase(ind_vars);
stdList = stdBase(ind_vars);

time_part1 = datetime(1992,10,2):datetime(2014,7,1,12,0,0);
time_part2 = datetime(2014,7,1, 12, 0, 0):hours(3):dateshift(datetime(datevec(now-1)), 'start', 'day');
time_pool = [time_part1(:); time_part2(:)];

if timeTick<time_pool(1) ||  timeTick>time_pool(end)
    error(['No available HYCOM data before 1992-10-02 or after ',datestr(now-1, 'yyyy-mm-dd'),'!'])
end

ind_rtime = wisefind(time_pool, timeTick);
time_hycom = time_pool(ind_rtime);

reg_names = abs(round(region));
geo_tag = ['W',num2str(reg_names(1)),'E',num2str(reg_names(2)),'S',num2str(reg_names(3)),'N',num2str(reg_names(4))];
aimfile = fullfile(aimpath, [geo_tag, '_',datestr(time_hycom,'yyyymmddTHHMMZ'),'.mat']);

%% Download
if exist(aimfile, 'file')~=0
    disp('It has been downloaded before')
    D = load(aimfile);
else
    if nargin< 5
        URL = get_URL(time_hycom);
    end
    %     ncdisp(URL)  % debug
    nc_dims = {'lon','lat','depth','time'};
    lonAll = ncread(URL, nc_dims{1});
    latAll = ncread(URL, nc_dims{2});
    depAll = ncread(URL, nc_dims{3});
    timeAll = datetime(datevec(ncread(URL, nc_dims{4})/24+datenum(2000,1,1)));

    if time_hycom<min(timeAll)-days(1) || time_hycom>max(timeAll)+days(1)
        error('The given time is outside the time range of the specified HYCOM product!')
    end
    if  min(lonAll)<0
        if max(region(1:2))>180
            ind_west = region>180;
            region(ind_west) = region(ind_west)-360;
        end
    end
    indLons = wisefind(lonAll, region(1:2));
    indLats = wisefind(latAll, region(3:4));
    indTime = wisefind(timeAll, time_hycom);

    D.lon = lonAll(min(indLons):max(indLons));
    D.lat = latAll(min(indLats):max(indLats));
    D.depth = depAll;
    D.time = timeAll(indTime);

    tole_time = days(1); % you can increase this value to replace missing data with data on adjacent days
    dev_time = D.time-timeTick;
    if abs(dev_time) > tole_time
        warning on
        warning(['HYCOM data is missing on ',datestr(time_hycom, 'yyyy-mm-dd') ,', and this function has stopped'])
        return
    end
    D.dev_time = dev_time;

    nVars = numel(stdList);
    nLayers = numel(D.depth);
    for iVar = 1:nVars
        stdName = stdList{iVar};
        varName = varList{iVar};
        if strcmp(varName, 'ssh')
            varData = squeeze(ncread(URL, stdName, [min(indLons),min(indLats),indTime], [abs(diff(indLons))+1,abs(diff(indLats))+1,1])); %#ok<*NASGU>
        else
            varData = squeeze(ncread(URL, stdName, [min(indLons),min(indLats),1,indTime], [abs(diff(indLons))+1,abs(diff(indLats))+1,nLayers,1]));
        end
        D.(varName) = varData;
        clear varData
    end
    save(aimfile, '-struct', 'D')
end

cst = toc;
disp(['It takes ', num2str(cst,'%.2f'),' secs to download ', datestr(time_hycom,'yyyymmddTHHMMZ'), '.mat'])
end

function URL = get_URL(timeTick)
% Automatically select a HYCOM product according to the given time
% DO NOT arbitrarily change the order below

if timeTick < datetime(1992,10,2)
    error('No available HYCOM data set before 1992-10-02!')

    % ------GLBv0.08 (2014-7-1 to 2020-2-19, 3-hourly, 40 levels, 0.08*0.08) --- High-priority
elseif timeTick >= datetime(2018,1,1,12,0,0) && timeTick <= datetime(2020,2,19,9,0,0)  % checked  0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0?';
elseif timeTick >= datetime(2017,10,1,12,0,0) && timeTick <= datetime(2018,3,20,9,0,0)  % checked  0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9?';
elseif timeTick >= datetime(2017,6,1,12,0,0) && timeTick <= datetime(2017,10,1,9,0,0)    % checked   -180-180, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7?';
elseif timeTick >= datetime(2017,2,1,12,0,0) && timeTick <= datetime(2017,6,1,9,0,0)      % checked  0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.8?';
elseif timeTick >= datetime(2016,5,1,12,0,0) && timeTick <= datetime(2017,2,1,9,0,0)       % checked   -180-180, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.2?';
elseif timeTick >= datetime(2014,7,1,12,0,0) && timeTick <= datetime(2016,9,30,9,0,0)     % checked   -180-180, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3?';

    % ------GLBy0.08 (2018-12-4 to ongoing, 3-hourly, 40 levels, 0.08*0.04)
elseif timeTick >= datetime(2018,12,4,12,0,0)                                                    % checked  0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?';

    % ------GLBu0.08 (1992-10-2 to 2018-11-20, daily, 40 levels, 0.08*0.08)
elseif timeTick >= datetime(2016,4,18) && timeTick <= datetime(2018,11,20)  % checked   0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2?';
elseif timeTick >= datetime(2014,4,7) && timeTick <= datetime(2016,4,18)      % checked    0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1?';
elseif timeTick >= datetime(2013,8,17) && timeTick <= datetime(2014,4,8)      % checked    0-360, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0?';
elseif timeTick >= datetime(1995,8,1) && timeTick <= datetime(2012,12,31)     % checked   -180-180, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1?';
elseif timeTick >= datetime(2012,1,25) && timeTick <= datetime(2013,8,20)    % checked    0-360, -80-80 (20120903-20121202 are missing)
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9?';
elseif timeTick >= datetime(1992,10,2) && timeTick <= datetime(1995,7,31)     % checked   -180-180, -80-80
    URL = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.0?';
end
expName = URL(end-18:end-1);
disp(['HYCOM_',expName, ' is being downloaded, please wait...'])
end

function indMin = wisefind(varBase, varFind)
% Find the closest index

varBase = varBase(:)';
varFind = varFind(:);
diff_vals = abs(varBase-varFind);  % This line might give an error when running on some older versions of MATLAB? I am not sure.
[~, indMin] = sort(diff_vals, 2);
indMin = indMin(:, 1);
end

%% All HYCOM products (global)
% GLBy0.08
% expt_93.0
%
% GLBv0.08
%  expt_93.0; expt_92.9; expt_57.7; expt_92.8; expt_57.2; expt_56.3;
%  expt_53.X (GOFS 3.1-Reanalysis) --------- this data set is not used.
%
% GLBa0.08 --------- this product has quite different coordinate system, and it is not used.
% expt_91.2; expt_91.1; expt_91.0; expt_90.9; expt_90.8; expt_90.6
%
% GLBu0.08
% expt_91.2; expt_91.1;expt_91.0; expt_90.9;
% expt_19.1; expt_19.0  (GOFS 3.0-Reanalysis)
%
%% Debug
% timeList = [datetime(1993,1,1) datetime(1996,1,1) datetime(2012, 3,1) datetime(2014, 1,1) ...
%     datetime(2014, 5,1) datetime(2015, 1,1) datetime(2016, 1,1) datetime(2017, 5,1) ...
%     datetime(2017, 9,1) datetime(2018, 2,1) datetime(2019, 3,1) datetime(2022, 5,1)];

