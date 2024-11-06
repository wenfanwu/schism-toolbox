function tripcolor(tri, x, y, z, varargin)
% Visualize the variable on triangular mesh
%
%% Syntax
% tripcolor(tri, x, y, z, varargin)
% tripcolor(tri, x, y, z)
%
%% Description
% tripcolor(tri, x, y, z) visualizes the variables defined at the
% triangular nodes or cells. 
% tripcolor(tri, x, y, z, varargin) specifies your own options.
%
%% Example
% 
%% Input Arguments
% tri --- the connectivity table of triangular mesh
% x --- x-direction coordinates
% y --- y-direction coordinates
% z --- variable defined at the nodes or cells; x, y, and z must have the
% same dimension.
% varagin --- options, the same as the 'patch' function
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: patch

%% Parse inputs
i34 = ~isnan(tri(:,end))+3;

if numel(find(i34==4)) == 0
    tri = tri(:, 1:3);
end

if size(tri,1) == length(z(:))
    opts_def = {'FaceColor','flat', 'EdgeAlpha', 0.1, 'CDataMapping','scaled', 'EdgeColor', 'none', 'LineWidth', 0.0025};
else
    opts_def = {'FaceColor','interp', 'EdgeAlpha', 0.1, 'CDataMapping','scaled', 'EdgeColor', 'none', 'LineWidth', 0.0025};
end

varargin = [opts_def(:)', varargin(:)'];
%% Display
patch('Faces',tri,'Vertices',[x(:) y(:)], 'FaceVertexCData',z(:), varargin{:})

end