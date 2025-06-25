function hotstart = write_schism_hotstart(Mobj, InitCnd, start_time)
% Write hotstart.nc for SCHISM model.
%
%% Syntax
% hotstart = write_schism_hotstart(Mobj, InitCnd, start_time)
%
%% Description
% hotstart = write_schism_hotstart(Mobj, InitCnd, start_time) writes the
%       hotstart.nc file and return the datastruct containing all variables.
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% InitCnd - initial condition; datastruct
%       the datastruct used to store initial fields.
% start_time - start time; datetime
%       the start time in hotstart.nc. Default: start_time = Mobj.time(1);
%
%% Output Arguments
% hotstart - hotstart data; datastruct
%       the datastruct containing all initial variables in hotstart.nc
%
%% Notes
% This function is under development to support additional tracer modules.
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 23 Jun 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_gr3

%% Parse inputs
if nargin < 3; start_time = Mobj.time(1); end
nTracers = Mobj.nTracers; maxLev = Mobj.maxLev; 
nElems = Mobj.nElems;  nNodes = Mobj.nNodes;  nEdges = Mobj.nEdges;

%% Reverse the vertical dimention
varList = {InitCnd.Variable};
for iVar = 1:numel(varList)
    varName = varList{iVar};
    ind_var = strcmp(varList, varName);
    varData = squeeze(InitCnd(ind_var).Data);
    
    % Ensure T/S values are valid
    switch varName
        case 'temp'; varData = max(-2, min(varData, 40));
        case 'salt'; varData = max(0, min(varData, 42));
    end
    % Last row represents the surface.
    if size(varData,1)==Mobj.maxLev  % match the vertical dimention
        InitCnd(ind_var).Data = flip(varData, 1);  % the first dimension is "depth"
    end
end
%% Extract data for all tracers
% Essential variables
hotstart.time = seconds(start_time-Mobj.time(1));   % start time
hotstart.eta2 = InitCnd(strcmp({InitCnd.Variable},'ssh')).Data;     % elevation
hotstart.tr_nd(1,:,:) = InitCnd(strcmp({InitCnd.Variable},'temp')).Data;  % temp
hotstart.tr_nd(2,:,:) = InitCnd(strcmp({InitCnd.Variable},'salt')).Data;   % salt

ind_mods = find(Mobj.tracer_counts~=0);
ind_mods(1:2) = [];  % the index of tracer modules (without temp & salt)

idx = 2; % at least two tracers (temp & salt)
for ii = ind_mods(:)'
    tracer_list = Mobj.tracer_sheet(3:end, ii);
    ind_tracers = find(~cellfun(@isempty, tracer_list));
    for jj = ind_tracers(:)'
        idx = idx+1;
        tracer_name = lower(tracer_list{jj});
        ind_var = strcmp({InitCnd.Variable}, tracer_name);
        varTmp = InitCnd(ind_var).Data;
        if min(size(varTmp))==1; varTmp = repmat(varTmp, Mobj.maxLev, 1); end
        
        if any(varTmp<0)
            warning on; warning(['Negative values exist in the ', tracer_name, ' data'])
            varTmp(varTmp<0) = 0;
        end
        hotstart.tr_nd(idx,:,:) = varTmp;
    end
end

% Variables defined at element centers
hotstart.tr_nd0 = hotstart.tr_nd;
for jj = 1:nTracers
    for kk = 1:maxLev
        hotstart.tr_el(jj,kk,:) = convert_schism_var(Mobj, squeeze(hotstart.tr_nd(jj, kk, :)), 'node2elem');
    end
end

% Misc variables (default)
hotstart.iths = 0;
hotstart.ifile = 1;
hotstart.idry_e = zeros(nElems,1);             % wet_dry flag at elements
hotstart.idry_s = zeros(nEdges,1);               % wet_dry flag at sides
hotstart.idry = zeros(nNodes,1);                % wet_dry flag
hotstart.we = zeros(maxLev, nElems);           % vertical velocity at elems
hotstart.su2 = zeros(maxLev, nEdges);           % side velocity u
hotstart.sv2 = zeros(maxLev, nEdges);           % side velocity v
hotstart.q2 = zeros(maxLev, nNodes);        
hotstart.xl = zeros(maxLev, nNodes);
hotstart.dfv = zeros(maxLev, nNodes);         % vertical diffusive coefficient
hotstart.dfh = zeros(maxLev, nNodes);         % horizontal diffusive coefficient
hotstart.dfq1 = zeros(maxLev, nNodes);
hotstart.dfq2 = zeros(maxLev, nNodes);
hotstart.nsteps_from_cold = 0;
hotstart.cumsum_eta = zeros(nNodes,1);

%% Begin to write (Hydro)
filepath = fullfile(Mobj.aimpath, 'hotstart.nc');
if exist(filepath,'file')==2; delete(filepath); end

nccreate(filepath,'time','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', hotstart.time);
nccreate(filepath,'iths','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'iths', hotstart.iths);
nccreate(filepath,'ifile','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'ifile', hotstart.ifile);

nccreate(filepath,'nsteps_from_cold','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')  % added since v5.9.0
ncwrite(filepath,'nsteps_from_cold', hotstart.nsteps_from_cold);
nccreate(filepath,'cumsum_eta','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'cumsum_eta', hotstart.cumsum_eta);

nccreate(filepath,'idry_e','Dimensions',{'elem', nElems},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry_e', hotstart.idry_e);
nccreate(filepath,'idry_s','Dimensions',{'side', nEdges},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry_s', hotstart.idry_s);
nccreate(filepath,'idry','Dimensions',{'node', nNodes},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry', hotstart.idry);

nccreate(filepath,'eta2','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'eta2', hotstart.eta2);
nccreate(filepath,'we','Dimensions',{'nVert', maxLev, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'we', hotstart.we);
nccreate(filepath,'tr_el','Dimensions',{'ntracers',nTracers, 'nVert', maxLev, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_el', hotstart.tr_el);
nccreate(filepath,'su2','Dimensions',{'nVert', maxLev, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'su2', hotstart.su2);
nccreate(filepath,'sv2','Dimensions',{'nVert', maxLev, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'sv2', hotstart.sv2);
nccreate(filepath,'tr_nd','Dimensions',{'ntracers',nTracers, 'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_nd', hotstart.tr_nd);
nccreate(filepath,'tr_nd0','Dimensions',{'ntracers',nTracers, 'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_nd0', hotstart.tr_nd0);
nccreate(filepath,'q2','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'q2', hotstart.q2);
nccreate(filepath,'xl','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'xl', hotstart.xl);
nccreate(filepath,'dfv','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfv', hotstart.dfv);
nccreate(filepath,'dfh','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfh', hotstart.dfh);
nccreate(filepath,'dfq1','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfq1', hotstart.dfq1);
nccreate(filepath,'dfq2','Dimensions',{'nVert', maxLev, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfq2', hotstart.dfq2);

%% USE_ICE
if strcmpi(Mobj.use_ice, 'yes')
    ntr_ice = Mobj.ntr_ice; % # of ice tracers (in order: 1: ice mass; 2: ice conc; 3: snow mass)
    
    hotstart.ice_free_flag = 1;     % start from ice-free conditions
    hotstart.ice_free_flag2 = 1;
    hotstart.ice_surface_T = -1.8*ones(nNodes, 1);
    hotstart.ice_water_flux = zeros(nNodes, 1);
    hotstart.ice_heat_flux = zeros(nNodes, 1);
    hotstart.ice_velocity_x = zeros(nNodes, 1);
    hotstart.ice_velocity_y = zeros(nNodes, 1);
    hotstart.ice_sigma11 = zeros(nElems, 1);
    hotstart.ice_sigma12 = zeros(nElems, 1);
    hotstart.ice_sigma22 = zeros(nElems, 1);
    hotstart.ice_ocean_stress = zeros(2, nNodes);
    hotstart.ice_tracers = ones(ntr_ice, nNodes);

    nccreate(filepath,'ice_free_flag','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_free_flag', hotstart.ice_free_flag);
    nccreate(filepath,'ice_free_flag2','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_free_flag2', hotstart.ice_free_flag2);

    ICE_vars = {'ice_free_flag', 'ice_free_flag2','ice_surface_T','ice_water_flux','ice_heat_flux','ice_velocity_x','ice_velocity_y'};
    for iVar = 3:numel(ICE_vars)
        varName = strtrim(ICE_vars{iVar});
        varData = hotstart.(varName);

        nccreate(filepath, varName,'Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
    
    ICE_vars2 = {'ice_sigma11','ice_sigma12','ice_sigma22'};
    for iVar = 1:numel(ICE_vars2)
        varName = strtrim(ICE_vars2{iVar});
        varData = hotstart.(varName);

        nccreate(filepath, varName,'Dimensions',{'elem', nElems},'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
    
    nccreate(filepath,'ice_ocean_stress','Dimensions',{'two',2, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_ocean_stress', hotstart.ice_ocean_stress);
    
    nccreate(filepath,'ice_tracers','Dimensions',{'ice_ntr', ntr_ice, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_tracers', hotstart.ice_tracers);
end
%% USE_COSINE
if strcmpi(Mobj.use_cosine, 'yes')
    ndelay = Mobj.ndelay;  % ndelay in the cosine.nml

    ind_cos = [7 10 8 9];  % bug
    hotstart.COS_sS2 = squeeze(hotstart.tr_el(ind_cos(1), :,:));
    hotstart.COS_sDN = squeeze(hotstart.tr_el(ind_cos(2), :,:));
    hotstart.COS_sZ1 = squeeze(hotstart.tr_el(ind_cos(3), :,:));
    hotstart.COS_sZ2 = squeeze(hotstart.tr_el(ind_cos(4), :,:));

    hotstart.COS_nstep = ones(maxLev, nElems);
    hotstart.COS_mS2 = permute(repmat(hotstart.COS_sS2, [1 1 ndelay]), [3 1 2]);
    hotstart.COS_mDN = permute(repmat(hotstart.COS_sDN, [1 1 ndelay]), [3 1 2]);
    hotstart.COS_mZ1 = permute(repmat(hotstart.COS_sZ1, [1 1 ndelay]), [3 1 2]);
    hotstart.COS_mZ2 = permute(repmat(hotstart.COS_sZ2, [1 1 ndelay]), [3 1 2]);

    COS_vars1 = {'COS_mS2', 'COS_mDN', 'COS_mZ1', 'COS_mZ2'};
    COS_vars2 = {'COS_sS2', 'COS_sDN', 'COS_sZ1', 'COS_sZ2'};
    for iVar = 1:numel(COS_vars1)
        varName = strtrim(COS_vars1{iVar});
        varData = hotstart.(varName);
        varData(varData<0) = 0;
        nccreate(filepath, varName,'Dimensions',{'seven', 7, 'nVert', maxLev, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);

        varName = strtrim(COS_vars2{iVar});
        varData = hotstart.(varName);
        varData(varData<0) = 0;
        nccreate(filepath, varName,'Dimensions',{'nVert', maxLev, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
     nccreate(filepath, 'COS_nstep','Dimensions',{'nVert', maxLev, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
     ncwrite(filepath, 'COS_nstep', hotstart.COS_nstep);
end

%% USE_SED3D
if strcmpi(Mobj.use_sed3d, 'yes')
    sed_class = Mobj.sed_class;  % # of sediment types
    nbed = Mobj.nbed;  % # of bed layers
    MBEDP = Mobj.MBEDP;

    hotstart.SED3D_dp = InitCnd.SED3D_dp;
    hotstart.SED3D_rough = InitCnd.SED3D_rough;
    hotstart.SED3D_bed = InitCnd.SED3D_bed;
    hotstart.SED3D_bedfrac = InitCnd.SED3D_bedfrac;

    nccreate(filepath,'SED3D_dp','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'SED3D_dp', hotstart.SED3D_dp);
    nccreate(filepath,'SED3D_rough','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'SED3D_rough', hotstart.SED3D_rough);
    nccreate(filepath,'SED3D_bed','Dimensions',{'MBEDP', MBEDP, 'nbed', nbed, 'elem', nElems},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'SED3D_bed', hotstart.SED3D_bed);
    nccreate(filepath,'SED3D_bedfrac','Dimensions',{'sed_class', sed_class, 'nbed', nbed, 'elem', nElems},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'SED3D_bedfrac', hotstart.SED3D_bedfrac);
end

%% USE_SED2D


%% USE_WWM



disp('hotstart.nc has been created successfully!')
end

% Notes:
% Initialize some arrays and constants
% tempmin=-2.d0; tempmax=40.d0; saltmin=0.d0; saltmax=42.d0
% pr1=0.d0; pr2=0.d0; pr=prmsl_ref !uniform pressure (the const. is unimportant)
% uthnd=-99.d0; vthnd=-99.d0; eta_mean=-99.d0; !uth=-99.d0; vth=-99.d0; !flags
% elevmax=-1.d34; dav_maxmag=-1.d0; dav_max=0.d0
% tr_el=0.d0
% timer_ns=0.d0
% iwsett=0; wsett=0.d0 !settling vel.
% rough_p=1.d-4 !>0 for SED
% tau_bot_node=0.d0
% uu2 = 0.0_rkind
% vv2 = 0.0_rkind
% ww2 = 0.0_rkind
% tr_nd0=0.d0
% istack0_schout=0
% cumsum_eta=0.d0
% nsteps_from_cold=0












