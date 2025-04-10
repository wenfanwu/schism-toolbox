function varBtm = get_schout_btm(Mobj, varTri, N)
% Get the variable at the near-bottom layer
% 
%% Syntax
% varBtm = get_schout_btm(Mobj, varTri, N)
%
%% Description 
% varBtm = get_schout_btm(Mobj, varTri, N) gets the near-bottom data
%
%% Examples
% varBtm = get_schout_btm(Mobj, varTri, 1)
%
%% Input Arguments
% Mobj --- the mesh object (with vertical layers added)
% N --- the N-th level from the bottom. N = 0 is the bottom layer.
%
%% Output Arguments
% varBtm --- variables at the N-th layer from the bottom.
% 
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2023-11-27.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
if nargin < 3
    N = 1;
end

ind_btm = sum(~isnan(Mobj.depLayers));
ind_btm = max(ind_btm-N, 1);
varBtm = arrayfun(@(x) varTri(ind_btm(x), x), 1:Mobj.nNodes);

end