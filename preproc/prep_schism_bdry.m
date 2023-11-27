function D = prep_schism_bdry(Mobj, dst)
% Prepare boundary inputs (Hydrological part)
%
%% Syntax
%
%
%% Description
%
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
% Created by Wenfan Wu, Ocean Univ. of China in 2022.
% Last Updated on 2022-05-17.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also:

%% Parse inputs
switch dst
    case 'hycom_bys'
        D = get_hycom_bdry_bys(Mobj); % Make sure the datapath inside it is valid.
    case 'hycom_clim'
        D = get_hycom_bdry_clim(Mobj);
end

end