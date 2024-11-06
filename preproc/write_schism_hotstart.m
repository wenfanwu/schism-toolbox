function Hs = write_schism_hotstart(Mobj, InitCnd, start_time)
% Write hotstart.nc for SCHISM (NOT Completed Yet)
%
%% Syntax
% hotstart_data = write_schism_hotstart(Mobj, InitCnd, start_time)
%
%% Description
% hotstart_data = write_schism_hotstart(Mobj, InitCnd, start_time) writes the
% hotstart.nc file and return a datastruct containing all the variables. 
%
%% Input Arguments
% Mobj --- the mesh object
% InitCnd --- the datastruct contains initial states.
% start_time --- the start time in hotstart.nc
%
%% Output Arguments
% hotstart_data --- a datastruct containing all the variables ready for the
% hotstart
%
%% Notes
% this function still needs a lot of improvements, as many of the initial
% variables are given initial values that may be not appropriate. 
% 
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2023-11-26. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: ncwrite

%% Parse inputs
nTracers = Mobj.nTracers;
nLevs = Mobj.maxLev;
nElems = Mobj.nElems;
nNodes = Mobj.nNodes;
nEdges = Mobj.nEdges;

%% The last row represents the surface in SCHISM
InitCnd.salt = flip(InitCnd.salt, 2);
InitCnd.temp = flip(InitCnd.temp, 2);

%% Make sure the T/S values are valid
temp_min = -2; temp_max = 40;
InitCnd.temp = max(temp_min, InitCnd.temp);
InitCnd.temp = min(temp_max, InitCnd.temp);

salt_min = 0; salt_max = 42;
InitCnd.salt = max(salt_min, InitCnd.salt);
InitCnd.salt = min(salt_max, InitCnd.salt);

%% NetCDF
Hs.time = seconds(start_time-Mobj.time(1));   % start time
Hs.eta2 = InitCnd.ssh;                                      % elevation
Hs.tr_nd(1,:,:) = InitCnd.temp';
Hs.tr_nd(2,:,:) = InitCnd.salt';

ind_mods = find(Mobj.tracer_counts~=0);
ind_mods(1:2) = [];  % index of the tracer model

sn = 2; % at least two tracers in the model (temp&salt)
for ii = ind_mods(:)'
    tracer_list = Mobj.tracer_sheet(3:end, ii);
    ind_tracers = find(~cellfun(@isempty, tracer_list));
    for jj = ind_tracers(:)'
        sn = sn+1;
        tracer_name = lower(tracer_list{jj});
        varTmp = InitCnd.(tracer_name);
        if min(size(varTmp))==1
            varTmp = repmat(varTmp, 1, Mobj.maxLev);
        end
        if sum(varTmp<0)~=0
            warning on
            warning(['Negative values exist in the ', tracer_name, ' data'])
            varTmp(varTmp<0) = 0;  % need changes
        end
        Hs.tr_nd(sn,:,:) = varTmp';
    end
end

Hs.tr_nd0 = Hs.tr_nd;
for jj = 1:nTracers
    for kk = 1:nLevs
        var_nd = squeeze(Hs.tr_nd(jj, kk, :));
        Hs.tr_el(jj,kk,:) = convert_schism_var(Mobj, var_nd, 'node2elem');
    end
end

Hs.iths = 0;
Hs.ifile = 1;
Hs.idry_e = zeros(nElems,1);             % wet_dry flag at elements
Hs.idry_s = zeros(nEdges,1);               % wet_dry flag at sides
Hs.idry = zeros(nNodes,1);                % wet_dry flag
Hs.we = zeros(nLevs, nElems);           % vertical velocity at elems
Hs.su2 = zeros(nLevs, nEdges);           % side velocity u
Hs.sv2 = zeros(nLevs, nEdges);           % side velocity v
Hs.q2 = zeros(nLevs, nNodes);        
Hs.xl = zeros(nLevs, nNodes);
Hs.dfv = zeros(nLevs, nNodes);         % vertical diffusive coefficient
Hs.dfh = zeros(nLevs, nNodes);         % horizontal diffusive coefficient
Hs.dfq1 = zeros(nLevs, nNodes);
Hs.dfq2 = zeros(nLevs, nNodes);
Hs.nsteps_from_cold = 0;
Hs.cumsum_eta = zeros(Mobj.nNodes,1);

%% NetCDF
fileName = [Mobj.aimpath, 'hotstart.nc'];
if exist(fileName,'file')==2; delete(fileName); end

nccreate(fileName,'time','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'time', Hs.time);
nccreate(fileName,'iths','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(fileName,'iths', Hs.iths);
nccreate(fileName,'ifile','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(fileName,'ifile', Hs.ifile);

nccreate(fileName,'nsteps_from_cold','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')  % added since v5.9.0
ncwrite(fileName,'nsteps_from_cold', Hs.nsteps_from_cold);
nccreate(fileName,'cumsum_eta','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'cumsum_eta', Hs.cumsum_eta);

nccreate(fileName,'idry_e','Dimensions',{'elem', nElems},'Datatype','int32','Format','netcdf4')
ncwrite(fileName,'idry_e', Hs.idry_e);
nccreate(fileName,'idry_s','Dimensions',{'side', nEdges},'Datatype','int32','Format','netcdf4')
ncwrite(fileName,'idry_s', Hs.idry_s);
nccreate(fileName,'idry','Dimensions',{'node', nNodes},'Datatype','int32','Format','netcdf4')
ncwrite(fileName,'idry', Hs.idry);

nccreate(fileName,'eta2','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'eta2', Hs.eta2);
nccreate(fileName,'we','Dimensions',{'nVert', nLevs, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'we', Hs.we);
nccreate(fileName,'tr_el','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'tr_el', Hs.tr_el);
nccreate(fileName,'su2','Dimensions',{'nVert', nLevs, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'su2', Hs.su2);
nccreate(fileName,'sv2','Dimensions',{'nVert', nLevs, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'sv2', Hs.sv2);
nccreate(fileName,'tr_nd','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'tr_nd', Hs.tr_nd);
nccreate(fileName,'tr_nd0','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'tr_nd0', Hs.tr_nd0);
nccreate(fileName,'q2','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'q2', Hs.q2);
nccreate(fileName,'xl','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'xl', Hs.xl);
nccreate(fileName,'dfv','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'dfv', Hs.dfv);
nccreate(fileName,'dfh','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'dfh', Hs.dfh);
nccreate(fileName,'dfq1','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'dfq1', Hs.dfq1);
nccreate(fileName,'dfq2','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(fileName,'dfq2', Hs.dfq2);

%% USE_ICE
if strcmpi(Mobj.use_ice, 'yes')
    ntr_ice = 3; % # of ice tracers (in order: 1: ice mass; 2: ice conc; 3: snow mass)
     
    Hs.ice_free_flag = 1;     % start from ice-free conditions
    Hs.ice_free_flag2 = 1;
    Hs.ice_surface_T = -1.8*ones(nNodes, 1);
    Hs.ice_water_flux = zeros(nNodes, 1);
    Hs.ice_heat_flux = zeros(nNodes, 1);
    Hs.ice_velocity_x = zeros(nNodes, 1);
    Hs.ice_velocity_y = zeros(nNodes, 1);
    Hs.ice_sigma11 = zeros(nElems, 1);
    Hs.ice_sigma12 = zeros(nElems, 1);
    Hs.ice_sigma22 = zeros(nElems, 1);
    Hs.ice_ocean_stress = zeros(2, nNodes);
    Hs.ice_tracers = ones(ntr_ice, nNodes);

    nccreate(fileName,'ice_free_flag','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(fileName,'ice_free_flag', Hs.ice_free_flag);
    nccreate(fileName,'ice_free_flag2','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(fileName,'ice_free_flag2', Hs.ice_free_flag2);

    ICE_vars = {'ice_free_flag', 'ice_free_flag2','ice_surface_T','ice_water_flux','ice_heat_flux','ice_velocity_x','ice_velocity_y'};
    for iVar = 3:numel(ICE_vars)
        varName = strtrim(ICE_vars{iVar});
        varData = Hs.(varName);

        nccreate(fileName, varName,'Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
        ncwrite(fileName, varName, varData);
    end
    
    ICE_vars2 = {'ice_sigma11','ice_sigma12','ice_sigma22'};
    for iVar = 1:numel(ICE_vars2)
        varName = strtrim(ICE_vars2{iVar});
        varData = Hs.(varName);

        nccreate(fileName, varName,'Dimensions',{'elem', nElems},'Datatype','double','Format','netcdf4')
        ncwrite(fileName, varName, varData);
    end
    
    nccreate(fileName,'ice_ocean_stress','Dimensions',{'two',2, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(fileName,'ice_ocean_stress', Hs.ice_ocean_stress);
    
    nccreate(fileName,'ice_tracers','Dimensions',{'ice_ntr', ntr_ice, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(fileName,'ice_tracers', Hs.ice_tracers);
end
%% USE_COSINE
if strcmpi(Mobj.use_cosine, 'yes')
    Mobj.ndelay = 7;
    ind_cos = [7 10 8 9];  % bug
    Hs.COS_sS2 = squeeze(Hs.tr_el(ind_cos(1), :,:));
    Hs.COS_sDN = squeeze(Hs.tr_el(ind_cos(2), :,:));
    Hs.COS_sZ1 = squeeze(Hs.tr_el(ind_cos(3), :,:));
    Hs.COS_sZ2 = squeeze(Hs.tr_el(ind_cos(4), :,:));

    Hs.COS_nstep = ones(nLevs, nElems);
    Hs.COS_mS2 = permute(repmat(Hs.COS_sS2, [1 1 Mobj.ndelay]), [3 1 2]);
    Hs.COS_mDN = permute(repmat(Hs.COS_sDN, [1 1 Mobj.ndelay]), [3 1 2]);
    Hs.COS_mZ1 = permute(repmat(Hs.COS_sZ1, [1 1 Mobj.ndelay]), [3 1 2]);
    Hs.COS_mZ2 = permute(repmat(Hs.COS_sZ2, [1 1 Mobj.ndelay]), [3 1 2]);

    COS_vars1 = {'COS_mS2', 'COS_mDN', 'COS_mZ1', 'COS_mZ2'};
    COS_vars2 = {'COS_sS2', 'COS_sDN', 'COS_sZ1', 'COS_sZ2'};
    for iVar = 1:numel(COS_vars1)
        varName = strtrim(COS_vars1{iVar});
        varData = Hs.(varName);
        varData(varData<0) = 0;
        nccreate(fileName, varName,'Dimensions',{'seven', 7, 'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(fileName, varName, varData);

        varName = strtrim(COS_vars2{iVar});
        varData = Hs.(varName);
        varData(varData<0) = 0;
        nccreate(fileName, varName,'Dimensions',{'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(fileName, varName, varData);
    end
     nccreate(fileName, 'COS_nstep','Dimensions',{'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
     ncwrite(fileName, 'COS_nstep', Hs.COS_nstep);
end

%% USE_SED2D


%% USE_SED3D


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












