function DS = get_hycom_bdry_clim(Mobj) %#ok<*STOUT>
% Extract climatology boundary inputs from HYCOM data
% 
%% Syntax
% 
%
%% Description 
% 
%
%% Examples
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
% Last Updated on 2022-10-25.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
filepath = 'hycom_clim_vars_1995_2020.mat';  % can be changed
load(filepath) %#ok<LOAD>

varList = {'ssh', 'temperature', 'salinity', 'uvel', 'vvel'};
nickList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};
RD = hycom_clim;
bdry_time = unique(dateshift(Mobj.time, 'start', 'months'));
%% Extract
depRaw = abs(RD.depth);
lon_bnd = Mobj.lon(Mobj.obc_nodes_tot);
lat_bnd = Mobj.lat(Mobj.obc_nodes_tot);

clear DS
for iVar = 1:5
    D.ind = Mobj.obc_nodes_tot;
    D.lon = Mobj.lon(D.ind);
    D.lat = Mobj.lat(D.ind);
    D.depth = depRaw;
    D.time = bdry_time;
    D.unit_time = 'month';

    varName = varList{iVar};
    nickName = nickList{iVar};
    varClim = RD.(varName);
    varClim = permute(varClim, [2 1 3 4]);  % lon*lat*depth*time
    varBnds = get_bnd_vars(lon_bnd, lat_bnd, RD.lon, RD.lat, varClim);  % nps*depth*time
    D.var = clim_full(varBnds, D.time);

    DS.(nickName) = D;  % nps*depth*time
    clear D
end
end

function varNew = clim_full(varRaw, timeRaw)
% Expand the variable matrix on its time dimension
[n1, n2, ~] = size(varRaw);

time_months = month(timeRaw);
month_ticks = sort(unique(time_months));
nt = numel(month_ticks);
varNew = nan(n1, n2, nt);

for iMon = month_ticks(:)'
    ind_mons = time_months==iMon;
    nMonths = numel(find(ind_mons));
    if nMonths~=0
        varNew(:,:,ind_mons) = repmat(varRaw(:,:,iMon), [1 1 nMonths]);
    end
end
end












