function time = any2time(timeRaw, time_interval, base_date)
% Convert various time formats to datetime
%
%% Syntax
%   time = any2time(timeRaw)
%   time = any2time(timeRaw, time_interval)
%   time = any2time(timeRaw, time_interval, base_date)
%
%% Description
%   time = any2time(timeRaw) converts various time formats to datetime
%   time = any2time(timeRaw, time_interval) specifies the time interval
%   time = any2time(timeRaw, time_interval, base_date) specifies the base date
% 
%% Example
%   timeRaw = datenum(2000:2010,1,1) - datenum(1900, 1, 1);
%   time = any2time(timeRaw, 'days', [1900 1 1]);
%
%% Input Arguments
%   timeRaw       - Time data (datenum, datestr, datevec, datetime)
%   time_interval - Time unit: 'days', 'hours', 'minutes', 'seconds' (optional)
%   base_date     - Reference date for numerical conversions (optional)
%
%% Output Arguments
%   time          - Time vector/matrix in datetime format.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 21 Feb 2025. 
% Email: wwu@vims.edu
%
% See also: datetime, datevec, datenum, datestr

%% Set Default Values
if nargin < 2, time_interval = 'days'; end
if nargin < 3, base_date = [0 0 0]; end

% Ensure base_date is in full datevec format
base_date_full = [1900 1 1 0 0 0]; 
base_date_full(1:length(base_date)) = base_date(:)';

if ischar(timeRaw) && isscalar(timeRaw)
    timeRaw = {timeRaw};
end

% Preserve input size for reshaping later
original_size = size(timeRaw);
if size(timeRaw, 2) == 6 && all(timeRaw(:,2) >= 1 & timeRaw(:,2) <= 12 & mod(timeRaw(:,2), 1) == 0)
    original_size(2) = 1;
else
    timeRaw = timeRaw(:);
end
%% Convert time based on input type
if isnumeric(timeRaw) % If datenum or datevec
    if size(timeRaw,2) == 6  % datevec
        time = datetime(timeRaw);
    else %  datenum
        switch lower(time_interval(1:2)) % Match first two letters
            case 'da', time_scale = 1;
            case 'ho', time_scale = 1/24;
            case 'mi', time_scale = 1/(24*60);
            case 'se', time_scale = 1/(24*60*60);
            otherwise, error('Invalid time interval: %s', time_interval);
        end
        time = datetime(datevec(timeRaw*time_scale + datenum(base_date_full))); %#ok<DATNM>
    end
elseif iscell(timeRaw) || isstring(timeRaw)% If datestr
    time = datetime(timeRaw);
elseif isdatetime(timeRaw) % If already datetime
    time = timeRaw;
else
    error('Unsupported time format: Input must be numeric, char, cell, string, or datetime.');
end

% Restore original shape
time = reshape(time, original_size);
end
