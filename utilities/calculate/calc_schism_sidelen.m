function [node_lens, elem_lens, edge_lens] = calc_schism_sidelen(Mobj)
% Calculate the side lengths (meters) of triangular cells
% 
%% Syntax
% [node_lens, elem_lens, edge_lens] = calc_schism_sidelen(Mobj)
%
%% Description 
% [node_lens, elem_lens, edge_lens] = calc_schism_sidelen(Mobj) calculates
% the side lengths of each triangular cell
%
%% Example
% [node_lens, elem_lens, edge_lens] = calc_schism_sidelen(Mobj)
%
%% Input Arguments
% Mobj --- the mesh object
%
%% Output Arguments
% node_lens --- side lengths@nodes
% elem_lens --- side lengths@elements
% edge_lens --- side lengths@edges
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-05-18.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
xxx = Mobj.lon(Mobj.tri);
yyy = Mobj.lat(Mobj.tri);
xxx2 = circshift(xxx, 1, 2);
yyy2 = circshift(yyy, 1, 2);

edge_lens = distance(yyy, xxx, yyy2, xxx2, [6378.137 0.0818191910428158])*1000;
elem_lens = mean(edge_lens,2);
node_lens = accumarray(Mobj.tri(:), edge_lens(:), [], @mean);

end