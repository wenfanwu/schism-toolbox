function disp_schism_var(Mobj, varData, varargin)
% Visualize the variable on triangular mesh
%
%% Syntax
% disp_schism_var(Mobj, varData)
% disp_schism_var(Mobj, varData, varargin)
%
%% Description
% disp_schism_var(Mobj, varData) visualizes the variable on a triangular mesh.
% disp_schism_var(Mobj, varData, varargin) modifies the mesh styles.
%
%% Examples 
% 
%% Input Arguments
% Mobj --- the mesh object
% varData --- variable vector defined at the nodes or elements.
% varargin --- this function shares the same options with the matlab
% built-in fuction 'patch.m'; some useful options include 'EdgeColor',
% 'EdgeAlpha', 'LineWidth', and so on. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 6 Jun. 2022. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: tripcolor and patch

%% Parse inputs
if numel(find(size(varData)~=1))~=1
    error('the input variable must be a vector!')
end
if numel(varData)~=Mobj.nNodes && numel(varData)~=Mobj.nElems
    error('the length of input variable is not consistent with the given mesh grid!')
end
varData = double(varData(:));

%% Display
tripcolor(Mobj.tri, Mobj.lon, Mobj.lat, varData, varargin{:})
colorbar
colormap(jet(25))
hold on; box on
xlabel('Longitude (°E)', 'FontWeight','bold')
ylabel('Latitude (°N)', 'FontWeight','bold')

end
