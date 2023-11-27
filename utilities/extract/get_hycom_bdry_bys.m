function DS = get_hycom_bdry_bys(Mobj) %#ok<*STOUT>
% HYCOM data in the Bohai, and Yellow Seas

%% Parse inputs
datapath = 'E:\Code-repository\Matlab-codes\functions-test\schism-toolbox-v1.0-beta\data\hycom\';  % NEED TO BE CHANGED

time_unit = 'days';
bdry_time = unique(dateshift(Mobj.time, 'start', time_unit));

cdate = datestr(bdry_time(1), 'yyyymmdd');
filepath = [datapath, cdate, '.mat'];
load(filepath, 'longitude', 'latitude','depth')
Hf.lon = longitude;
Hf.lat = latitude;
Hf.depth = abs(depth);
Hf.time = bdry_time;

nLons = numel(Hf.lon);
nLats = numel(Hf.lat);
nDeps = numel(Hf.depth);
nTimes = numel(Hf.time);
nNodes_obc = Mobj.nNodes_obc;

%% Interpolation
varList = {'ssh', 'temperature', 'salinity', 'uvel', 'vvel'}; % variable names in the original dataset
nickList = {'ssh', 'temp', 'salt', 'uvel', 'vvel'};  % standard variable names

nVars = numel(varList);
for iVar = 1:nVars

    D.ind = Mobj.obc_nodes_tot;
    D.lon = Mobj.lon(D.ind);
    D.lat = Mobj.lat(D.ind);
    D.depth = abs(Hf.depth);
    D.time = bdry_time;
    D.unit_time = 'day';

    indLon = minfind(Hf.lon, D.lon);
    indLat = minfind(Hf.lat, D.lat);
    indBnd = sub2ind([nLats, nLons], indLat, indLon);

    varName = varList{iVar};
    nickName = nickList{iVar};
    if strcmp(nickName, 'ssh')
        varAll = zeros(nNodes_obc, 1, nTimes);
    else
        varAll = zeros(nNodes_obc,nDeps, nTimes);
    end
    for iTime = 1:nTimes
        progressbar(iTime/nTimes)
        cdate = datestr(bdry_time(iTime), 'yyyymmdd');
        filepath = [datapath, cdate, '.mat'];
        if exist(filepath, 'file') ~=2  % for the missing files
            cdate2 = datestr(bdry_time(iTime)-days(1), 'yyyymmdd');
            filepath = [datapath, cdate2, '.mat'];
        end
        load(filepath, varName)
        varTmp = eval(varName);

        if strcmp(nickName, 'ssh')
            varAll(:, 1, iTime) = varTmp(indBnd);
        else
            varTmp = reshape(varTmp, [nLats*nLons nDeps]);
            varAll(:,:,iTime) = varTmp(indBnd, :);
        end
    end
    D.var = varAll;

    DS.(nickName) = D;
    clear D varAll
end

end

























































