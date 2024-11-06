function auto_center(base_size, adjust_factor)
% Automatically center the figure with minimal blank space around the plot
%
%% Syntax
% auto_center
% auto_center(base_size)
% auto_center(base_size, adjust_factor)
%
%% Description
% auto_center centerizes the current figure
% auto_center(base_size) specifies the base size of figure
% auto_center(base_size, adjust_factor) adjusts the figure height
% 
%% Examples 
% figure
% disp_schism_var(Mobj, Mobj.depth)
% auto_center
%
%% Input Arguments
% base_size - base size; double; 
%       specifies the base size of figure. base_size = 0.5 (default);
%       base_size >0 and <1;
% adjust_factor - adjust factor; double
%       adjusts the figure height slightly, adjust_factor = 1 (default);
%
%% Output Arguments
% None
% 
%% Notes
% This function was generated with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 29 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if nargin==0
    base_size = 0.5;
end
if nargin<2
    adjust_factor = 1;
end

%% Center figure
ax = gca;

xdar  = ax.DataAspectRatio(1);
ydar  = ax.DataAspectRatio(2);

dx = abs(diff(ax.XLim))/xdar;
dy = abs(diff(ax.YLim))/ydar;
aspect_ratio = dx/dy*ax.InnerPosition(4)/ax.OuterPosition(4);

fig_width=base_size*aspect_ratio;
fig_height=base_size*adjust_factor;

R = fig_height/fig_width;

if R>=1   % portrait
    fig_width = min(fig_width,1);
    fig_height = fig_width*R;
else  % landscape
    fig_height = min(fig_height,1);
    fig_width = fig_height/R;
end

left_margin = (1-fig_width)/2;
bottom_margin = (1-fig_height)/2;

set(gcf, 'Units', 'normalized', 'Position', [left_margin, bottom_margin, fig_width, fig_height]);

end
