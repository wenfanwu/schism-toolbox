function auto_center(target_ratio)
% Automatically center the figure
%
%% Syntax
% auto_center()
%
%% Description
% auto_center() centers the figure automatically
%
%% Examples 
% figure
% pcolor(peaks(20))
% shading flat
% auto_center
%
%% Input Arguments
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
% Last Updated on 2 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse input
if nargin == 0
    target_ratio = 0.618;  % the ratio of the visually "long" side to the screen
end
%% Center the figure
% get the size of screen (the main screen)
screen_size = get(0, 'ScreenSize');  % [left bottom width height]
screen_width_px = screen_size(3);
screen_height_px = screen_size(4);

% the # of pixels per inch of the screen (dpi)
dpi = get(0, 'ScreenPixelsPerInch');  % usually 96，but it may vary from PCs

% convert to inches (absolute measures)
screen_width_in = screen_width_px / dpi;
screen_height_in = screen_height_px / dpi;
R = screen_width_in/screen_height_in;

fig = gcf; ax = gca;

% set the units to normalized
set(fig, 'Units', 'normalized');
set(ax, 'Units', 'normalized');

% get the position of current figure and PlotBoxAspectRatio
pos = get(ax, 'Position');
pb_aspect = get(ax, 'PlotBoxAspectRatio');
pb_ratio = pb_aspect(2)/pb_aspect(1);

% calculates the visual aspect raiot (height/width)
visual_ratio = (pos(3)/pos(4)) * pb_ratio;

if pb_ratio >= 1
    % portrait
    fig_height = target_ratio;
    fig_width = (screen_height_in*fig_height/visual_ratio)/screen_width_in;
    fig_width = fig_width*1.1;
    
    if R<=1 && fig_width>0.5 % portrait screen
        fig_width = 0.5;
        fig_height = (screen_width_in*fig_width * visual_ratio)/screen_height_in;
    end
else
    % landscape
    fig_width = target_ratio;
    fig_height = (screen_width_in*fig_width * visual_ratio)/screen_height_in;
    fig_height = fig_height*1.1;

    if R>=1 && fig_height>0.5 % % portrait screen
        fig_height = 0.5;
        fig_width = (screen_height_in*fig_height/visual_ratio)/screen_width_in;
    end
end

% center the figure
left_margin = (1 - fig_width) / 2;
bottom_margin = (1 - fig_height) / 2;

set(fig, 'Position', [left_margin, bottom_margin, fig_width, fig_height], 'Color', 'w');

end
