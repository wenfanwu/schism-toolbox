function disp_schism_vgrid(Mobj, sect_info)
% Visualize the transect layers
%
%% Syntax
% disp_schism_vgrid(Mobj, SectNodes)
%
%% Description 
% disp_schism_vgrid(Mobj, SectNodes) visualizes the vertical layers on a
% given transect.
%
%% Example
%
%
%% Input Arguments
%
%
%% Output Arguments
% 
% 
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2022-05-18.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
depTrans = sect_info.depth;
nNodes = length(sect_info.lon);
depMax = min(depTrans, [], 1);

ind_btm = sum(~isnan(depTrans));
ind_stair = diff([ind_btm ind_btm(end)]);
for ii = 1:nNodes
    if ind_stair(ii) > 0
        depTrans(ind_btm(ii)+1:ind_btm(ii)+ind_stair(ii), ii) = depMax(ii);
    end
    if ind_stair(ii) < 0
        depTrans(ind_btm(ii)+ind_stair(ii)+1:ind_btm(ii), ii+1) = depMax(ii+1);
    end
end

%% Display
figure('Color', 'w')
disp_schism_hgrid(Mobj, [1 0])
colormap(jet(25))
hold on
scatter(sect_info.lon, sect_info.lat, 10, 'filled','k')
colorbar off
grid off; box off
axis off

figure('Color','w')
hold on
for iDep = 1:size(depTrans,1)
    plot(sect_info.dist, depTrans(iDep,:), 'Color','k')
    xlim([0 max(sect_info.dist)])
end
for iNode = 1:nNodes
    layNums = numel(find(~isnan(depTrans(:,iNode))));
    plot(sect_info.dist(iNode)*ones(1, layNums), depTrans(1:layNums,iNode), 'Color','k')
end
plot(sect_info.dist, depMax, 'Color','k', 'LineWidth',2)
xlabel('Along transect distance (km)');
ylabel('Depth (m)');

end










