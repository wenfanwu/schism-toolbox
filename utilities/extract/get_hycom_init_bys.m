function DS = get_hycom_init_bys(Mobj) 
% Get the hycom data as intial inputs
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
% Last Updated on 2022-05-20.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs 
init_time = Mobj.time(1);

datapath = 'E:\Code-repository\Matlab-codes\functions-test\schism-toolbox-v1.0-beta\data\hycom\';  % can be changed

cdate = datestr(init_time, 'yyyymmdd');
filepath = [datapath, cdate, '.mat'];

RD = load(filepath);

varList1 = {'ssh','salinity','temperature'};
varList2 = {'ssh','salt','temp'};
dim_order = [2 1 3];

nVars = numel(varList1);
for iVar = 1:nVars
    rawName = varList1{iVar};
    newName = varList2{iVar};
    D.lon = RD.longitude;
    D.lat = RD.latitude;
    D.depth = RD.depth;
    D.time = init_time;
    varRaw = squeeze(RD.(rawName));
    switch ndims(varRaw)
        case 2
            D.var = permute(varRaw, dim_order(1:2));
        case 3
            D.var = permute(varRaw, dim_order);
    end
    DS.(newName) = D;
    clear D
end

end
