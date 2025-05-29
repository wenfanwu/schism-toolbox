function auto_center(monitor_id, target_ratio)
% Automatically center the figure on the specified or current mouse screen
%
%% Syntax
% auto_center()
% auto_center(monitor_id)
% auto_center(monitor_id, target_ratio)
%
%% Description
% auto_center() centers the figure automatically.
% auto_center(monitor_id) specifies the monitor id.
% auto_center(monitor_id, target_ratio) specifies the aspect ratio.
%
%% Examples 
% figure
% pcolor(peaks(20))
% shading flat
% auto_center
%
%% Input Arguments
% monitor_id - the monitor id (optional); numeric
%       the monitor used to display the figure; If monitor_id is not
%       provided, the function automatically detects the monitor where the
%       mouse cursor is currently located.
% target_ratio - the ratio of long side (optional); numeric
%       the ratio of the visually longest side to the screen side.
%       Default: target_ratio = 0.618.
%
%% Output Arguments
% None
%
%% Notes
% This function works not so well for multiple-subplot cases.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2023. 
% Last Updated on 23 May 2025.
% Email: wwu@vims.edu
% 
% See also: 

%% Parse input
if nargin < 2 || isempty(target_ratio)
    target_ratio = 0.618;
end

monitors = get(0, 'MonitorPositions');  % Each row: [left bottom width height]
num_monitors = size(monitors, 1);

% Auto-detect monitor based on mouse position
if nargin < 1 || isempty(monitor_id)
    mouse_pos = get(0, 'PointerLocation');  % [x, y]
    monitor_id = find_monitor_by_position(mouse_pos, monitors);
end

if monitor_id < 1 || monitor_id > num_monitors
    error('Invalid monitor_id. You have %d monitors.', num_monitors);
end

target_monitor = monitors(monitor_id, :);  % [left bottom width height]
screen_width_px = target_monitor(3);
screen_height_px = target_monitor(4);

dpi = get(0, 'ScreenPixelsPerInch');
screen_width_in = screen_width_px / dpi;
screen_height_in = screen_height_px / dpi;
R = screen_width_in / screen_height_in;

fig = gcf; ax = gca;
set(fig, 'Units', 'normalized');
set(ax, 'Units', 'normalized');

pos = get(ax, 'Position');
pb_aspect = get(ax, 'PlotBoxAspectRatio');
pb_ratio = pb_aspect(2) / pb_aspect(1);
visual_ratio = (pos(3)/pos(4)) * pb_ratio;

if pb_ratio >= 1
    fig_height = target_ratio;
    fig_width = (screen_height_in*fig_height/visual_ratio)/screen_width_in;
    fig_width = fig_width*1.1;

    if R <= 1 && fig_width > 0.5
        fig_width = 0.5;
        fig_height = (screen_width_in*fig_width * visual_ratio)/screen_height_in;
    end
else
    fig_width = target_ratio;
    fig_height = (screen_width_in*fig_width * visual_ratio)/screen_height_in;
    fig_height = fig_height*1.1;

    if R >= 1 && fig_height > 0.5
        fig_height = 0.5;
        fig_width = (screen_height_in*fig_height/visual_ratio)/screen_width_in;
    end
end

left_margin = (1 - fig_width) / 2;
bottom_margin = (1 - fig_height) / 2;

% Convert to pixel coordinates
fig_left_px = target_monitor(1) + screen_width_px * left_margin;
fig_bottom_px = target_monitor(2) + screen_height_px * bottom_margin;
fig_width_px = screen_width_px * fig_width;
fig_height_px = screen_height_px * fig_height;

set(fig, 'Units', 'pixels');
set(fig, 'Position', [fig_left_px, fig_bottom_px, fig_width_px, fig_height_px], ...
         'Color', 'w');
end

function idx = find_monitor_by_position(pos, monitors)
% Find which monitor contains the position [x, y]
for i = 1:size(monitors, 1)
    m = monitors(i, :);
    if pos(1) >= m(1) && pos(1) <= m(1)+m(3) && ...
       pos(2) >= m(2) && pos(2) <= m(2)+m(4)
        idx = i;
        return;
    end
end
% If not found, fallback to monitor 1
idx = 1;
end
