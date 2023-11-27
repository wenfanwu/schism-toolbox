function auto_center(zoom_factor, aspect_ratio)
% Automatically center the figure
%
%% Syntax
% auto_center
% auto_center(zoom_factor)
% auto_center(zoom_factor, aspect_ratio)
%
%% Description
% auto_center() automatically centers the figure
% auto_center(zoom_factor) specifies the adjustment scale (>0)
% auto_center(zoom_factor, aspect_ratio) specifies aspect ratio (>0)
%
%% Example
% figure
% pcolor(peaks(20))
% auto_center
%
% figure
% subplot(211)
% pcolor(peaks(20))
% subplot(212)
% pcolor(peaks(20))
% auto_center(0.65, 0.5)
% 
%% Input Arguments
% zoom_factor --- zoom factor (>0); Default: zoom_factor = 0.5;
% aspect_ratio --- aspect ratio (>0), with higher value denoting longer x-axis
% visually; Default: aspect_ratio = 1; 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-04-14. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: get and set

%% Parse inputs
if nargin < 1
    zoom_factor = 0.5;
end
if nargin < 2
    aspect_ratio = 1;
end
if isempty(zoom_factor)
    zoom_factor = 0.5;
end

scnsize = get(0, 'ScreenSize');
H = scnsize(4);
L = scnsize(3);

dx = abs(diff(xlim));
dy = abs(diff(ylim));
%% Adjustment
udx = dx/L*aspect_ratio;
udy = dy/H;

pt = (1-zoom_factor)/2; 
if udx>=udy
    set(gcf, 'Units', 'normalized', 'Position', [pt, pt, zoom_factor, zoom_factor/udx*udy])
else
    set(gcf, 'Units', 'normalized', 'Position', [pt, pt, zoom_factor/udy*udx, zoom_factor])
end
end