function time = any2time(timeRaw, time_interval, base_date)
% Convert to datetime from any date formats
%
%% Syntax
% time = any2time(timeRaw)
% time = any2time(timeRaw, time_interval)
% time = any2time(timeRaw, time_interval, base_date)
%
%% Description
% time = any2time(timeRaw) converts the input time vector into the format
% of datetime.
%
% time = any2time(timeRaw, time_interval) specifies the time resolution of
% input time.
% time = any2time(timeRaw, time_interval, base_date) specifies the base date.
%
%% Example
% timeRaw = datenum(2000,1,1) - datenum(1900, 1, 1);
% time = any2time(timeRaw, 'days', [1900 1 1]);
% 
%% Input Arguments
% timeRaw --- input time vector with any date formats (datenum/datestr/datevec/datetime)
% 
% time_interval --- time interval,  it can be 'days', 'hours', 'minutes', or
% 'seconds'. Case is ignored, only the first two digits will be used.
% 
% base_date --- base date, e.g. base_date = [1900 1 1]; Default:
% base_date = [1900 1 1 0 0 0];
%
%% Output Arguments
% time --- time vector with datetime format.
%
%% Notes
% this function aims to convert a time vector of any date formats into
% 'datetime' format. Sometimes you may need to input the time
% origin, if your time vector with 'datenum' format starts at specific
% time origin. The time scale should be specified, if the provided vector is not
% daily.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2020. 
% Last Updated on 29 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: datetime

%% Parse inputs
switch nargin
    case 1
        time_interval = 'days';
        base_date = [0,0,0];
    case 2
        base_date = [0,0,0];
end
size_list = size(timeRaw);
timeRaw = timeRaw(:);

base_date_used = [1900 1 1 0 0 0];
base_date_used(1:length(base_date)) = base_date(:)';
%--------- From datenum
if isa(timeRaw,'numeric') && min(size(timeRaw,1),size(timeRaw,2)) == 1
    timeRaw = double(timeRaw);
    if strncmpi(time_interval,'months',2)
        time_gap = 1/31;
    end
    if strncmpi(time_interval,'days',2)
        time_gap = 1;
    end
    if strncmpi(time_interval,'hours',2)
        time_gap = 24;
    end
    if strncmpi(time_interval,'minutes',2)
        time_gap = 24*60;
    end
    if strncmpi(time_interval,'seconds',2)
        time_gap = 24*60*60;
    end
    time = datetime(datevec(timeRaw/time_gap+datenum(base_date_used)));
end
%--------- From datevec
if isa(timeRaw,'numeric') && (size(timeRaw,1) == 6 || size(timeRaw,2) == 6)
    timeRaw = double(timeRaw);
    time = datetime(timeRaw);
end
%--------- From datestr
if isa(timeRaw,'char')
    time = datetime(datevec(timeRaw));
end
%--------- From datetime
if isa(timeRaw,'datetime')
    time = timeRaw;
end

time = reshape(time, size_list);
end






