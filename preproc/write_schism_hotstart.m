function C = write_schism_hotstart(Mobj, InitCnd, start_time)
% Write hotstart.nc for SCHISM (Not completed yet)
%
%% Syntax
% C = write_schism_hotstart(Mobj, InitCnd, start_time)
%
%% Description
% C = write_schism_hotstart(Mobj, InitCnd, start_time) writes the
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
if nargin < 3; start_time = Mobj.time(1); end

nTracers = Mobj.nTracers; nLevs = Mobj.maxLev; 
nElems = Mobj.nElems; nNodes = Mobj.nNodes; nEdges = Mobj.nEdges;
%% Reverse 3-D variables
varList = {InitCnd.Variable};
for iVar = 1:numel(varList)
    varName = varList{iVar};
    ind_var = strcmp({InitCnd.Variable}, varName);
    varData = squeeze(InitCnd(ind_var).Data);
    
    % make sure the T/S values are valid
    switch varName
        case 'temp'; varData = max(-2, min(varData, 40));
        case 'salt'; varData = max(0, min(varData, 42));
    end
    % the last row represents the surface.
    if size(varData,1)==Mobj.maxLev
        InitCnd(ind_var).Data = flip(varData, 1);  % the first dimension is "depth"
    end
end
%% Extract data for all tracers
% essential variables
C.time = seconds(start_time-Mobj.time(1));   % start time
C.eta2 = InitCnd(strcmp({InitCnd.Variable},'ssh')).Data;     % elevation
C.tr_nd(1,:,:) = InitCnd(strcmp({InitCnd.Variable},'temp')).Data;  % temp
C.tr_nd(2,:,:) = InitCnd(strcmp({InitCnd.Variable},'salt')).Data;   % salt

ind_mods = find(Mobj.tracer_counts~=0);
ind_mods(1:2) = [];  % the index of tracer modules

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
        C.tr_nd(idx,:,:) = varTmp;
    end
end

% Variables defined at element centers
C.tr_nd0 = C.tr_nd;
for jj = 1:nTracers
    for kk = 1:nLevs
        C.tr_el(jj,kk,:) = convert_schism_var(Mobj, squeeze(C.tr_nd(jj, kk, :)), 'node2elem');
    end
end

% Assign variables
C.iths = 0;
C.ifile = 1;
C.idry_e = zeros(nElems,1);             % wet_dry flag at elements
C.idry_s = zeros(nEdges,1);               % wet_dry flag at sides
C.idry = zeros(nNodes,1);                % wet_dry flag
C.we = zeros(nLevs, nElems);           % vertical velocity at elems
C.su2 = zeros(nLevs, nEdges);           % side velocity u
C.sv2 = zeros(nLevs, nEdges);           % side velocity v
C.q2 = zeros(nLevs, nNodes);        
C.xl = zeros(nLevs, nNodes);
C.dfv = zeros(nLevs, nNodes);         % vertical diffusive coefficient
C.dfh = zeros(nLevs, nNodes);         % horizontal diffusive coefficient
C.dfq1 = zeros(nLevs, nNodes);
C.dfq2 = zeros(nLevs, nNodes);
C.nsteps_from_cold = 0;
C.cumsum_eta = zeros(Mobj.nNodes,1);

%% Begin to write (Hydro)
filepath = fullfile(Mobj.aimpath, 'hotstart.nc');
if exist(filepath,'file')==2; delete(filepath); end

nccreate(filepath,'time','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', C.time);
nccreate(filepath,'iths','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'iths', C.iths);
nccreate(filepath,'ifile','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'ifile', C.ifile);

nccreate(filepath,'nsteps_from_cold','Dimensions',{'one', 1},'Datatype','int32','Format','netcdf4')  % added since v5.9.0
ncwrite(filepath,'nsteps_from_cold', C.nsteps_from_cold);
nccreate(filepath,'cumsum_eta','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'cumsum_eta', C.cumsum_eta);

nccreate(filepath,'idry_e','Dimensions',{'elem', nElems},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry_e', C.idry_e);
nccreate(filepath,'idry_s','Dimensions',{'side', nEdges},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry_s', C.idry_s);
nccreate(filepath,'idry','Dimensions',{'node', nNodes},'Datatype','int32','Format','netcdf4')
ncwrite(filepath,'idry', C.idry);

nccreate(filepath,'eta2','Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'eta2', C.eta2);
nccreate(filepath,'we','Dimensions',{'nVert', nLevs, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'we', C.we);
nccreate(filepath,'tr_el','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'elem', nElems},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_el', C.tr_el);
nccreate(filepath,'su2','Dimensions',{'nVert', nLevs, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'su2', C.su2);
nccreate(filepath,'sv2','Dimensions',{'nVert', nLevs, 'side', nEdges},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'sv2', C.sv2);
nccreate(filepath,'tr_nd','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_nd', C.tr_nd);
nccreate(filepath,'tr_nd0','Dimensions',{'ntracers',nTracers, 'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'tr_nd0', C.tr_nd0);
nccreate(filepath,'q2','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'q2', C.q2);
nccreate(filepath,'xl','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'xl', C.xl);
nccreate(filepath,'dfv','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfv', C.dfv);
nccreate(filepath,'dfh','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfh', C.dfh);
nccreate(filepath,'dfq1','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfq1', C.dfq1);
nccreate(filepath,'dfq2','Dimensions',{'nVert', nLevs, 'node', nNodes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'dfq2', C.dfq2);

%% USE_ICE
if strcmpi(Mobj.use_ice, 'yes')
    ntr_ice = 3; % # of ice tracers (in order: 1: ice mass; 2: ice conc; 3: snow mass)
     
    C.ice_free_flag = 1;     % start from ice-free conditions
    C.ice_free_flag2 = 1;
    C.ice_surface_T = -1.8*ones(nNodes, 1);
    C.ice_water_flux = zeros(nNodes, 1);
    C.ice_heat_flux = zeros(nNodes, 1);
    C.ice_velocity_x = zeros(nNodes, 1);
    C.ice_velocity_y = zeros(nNodes, 1);
    C.ice_sigma11 = zeros(nElems, 1);
    C.ice_sigma12 = zeros(nElems, 1);
    C.ice_sigma22 = zeros(nElems, 1);
    C.ice_ocean_stress = zeros(2, nNodes);
    C.ice_tracers = ones(ntr_ice, nNodes);

    nccreate(filepath,'ice_free_flag','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_free_flag', C.ice_free_flag);
    nccreate(filepath,'ice_free_flag2','Dimensions',{'one', 1},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_free_flag2', C.ice_free_flag2);

    ICE_vars = {'ice_free_flag', 'ice_free_flag2','ice_surface_T','ice_water_flux','ice_heat_flux','ice_velocity_x','ice_velocity_y'};
    for iVar = 3:numel(ICE_vars)
        varName = strtrim(ICE_vars{iVar});
        varData = C.(varName);

        nccreate(filepath, varName,'Dimensions',{'node', nNodes},'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
    
    ICE_vars2 = {'ice_sigma11','ice_sigma12','ice_sigma22'};
    for iVar = 1:numel(ICE_vars2)
        varName = strtrim(ICE_vars2{iVar});
        varData = C.(varName);

        nccreate(filepath, varName,'Dimensions',{'elem', nElems},'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
    
    nccreate(filepath,'ice_ocean_stress','Dimensions',{'two',2, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_ocean_stress', C.ice_ocean_stress);
    
    nccreate(filepath,'ice_tracers','Dimensions',{'ice_ntr', ntr_ice, 'node', nNodes},'Datatype','double','Format','netcdf4')
    ncwrite(filepath,'ice_tracers', C.ice_tracers);
end
%% USE_COSINE
if strcmpi(Mobj.use_cosine, 'yes')
    Mobj.ndelay = 7;
    ind_cos = [7 10 8 9];  % bug
    C.COS_sS2 = squeeze(C.tr_el(ind_cos(1), :,:));
    C.COS_sDN = squeeze(C.tr_el(ind_cos(2), :,:));
    C.COS_sZ1 = squeeze(C.tr_el(ind_cos(3), :,:));
    C.COS_sZ2 = squeeze(C.tr_el(ind_cos(4), :,:));

    C.COS_nstep = ones(nLevs, nElems);
    C.COS_mS2 = permute(repmat(C.COS_sS2, [1 1 Mobj.ndelay]), [3 1 2]);
    C.COS_mDN = permute(repmat(C.COS_sDN, [1 1 Mobj.ndelay]), [3 1 2]);
    C.COS_mZ1 = permute(repmat(C.COS_sZ1, [1 1 Mobj.ndelay]), [3 1 2]);
    C.COS_mZ2 = permute(repmat(C.COS_sZ2, [1 1 Mobj.ndelay]), [3 1 2]);

    COS_vars1 = {'COS_mS2', 'COS_mDN', 'COS_mZ1', 'COS_mZ2'};
    COS_vars2 = {'COS_sS2', 'COS_sDN', 'COS_sZ1', 'COS_sZ2'};
    for iVar = 1:numel(COS_vars1)
        varName = strtrim(COS_vars1{iVar});
        varData = C.(varName);
        varData(varData<0) = 0;
        nccreate(filepath, varName,'Dimensions',{'seven', 7, 'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);

        varName = strtrim(COS_vars2{iVar});
        varData = C.(varName);
        varData(varData<0) = 0;
        nccreate(filepath, varName,'Dimensions',{'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
        ncwrite(filepath, varName, varData);
    end
     nccreate(filepath, 'COS_nstep','Dimensions',{'nVert', nLevs, 'elem', nElems}, 'Datatype','double','Format','netcdf4')
     ncwrite(filepath, 'COS_nstep', C.COS_nstep);
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












