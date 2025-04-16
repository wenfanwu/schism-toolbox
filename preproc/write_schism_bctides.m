function write_schism_bctides(Mobj, TideForc, bc_flags)
% Write the bctides.in file for SCHISM.
% 
%% Syntax
% write_schism_bctides(Mobj, TideForc, bc_flags)
%
%% Description 
% write_schism_bctides(Mobj, TideForc, bc_flags) creates bctide.in file
%
%% Examples
% tide_list = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
% TideForc = get_fes2014_tide(Mobj, tide_list);   
% write_schism_bctides(Mobj, TideForc, [5 5 4 4])
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct containing mesh info.
% TideForc - tidal forcing; datastruct
%       the datastruct containing tidal forcing data, which can be obtained
%       from the 'get_fes2014_tide.m' function. 
% bc_flags - boundary condition flags; double
%       the open boundary types specified for different tracers. bc_flags
%       is of M*N, where M is the # of open boundaries. M can be 1 in any
%       cases, meaning that the treatments of all tracers are identical on
%       all boundaries. In this way, you can specify different flags for
%       different boundaries. N indicates the # of activated tracer
%       modules, and N should be at least 4 (for purely hydrodynamic case).
%
% required fields in TideForc:
%       mandatory: tide_list; cutoff_depth (tide_list = {} means nTides=0, viz. no tidal forcing)
%       nTides~=0: nodal_factor; eq_arg;  
%
%       elev_flag=2: const_elev;  
%       elev_flag=[3,5]: elev_amp, elev_pha; 
%
%       uv_flag=2: const_flow; 
%       uv_flag=-4: nf_inflow; nf_outflow; 
%       uv_flag=[3,5]: u_amp; v_amp; u_pha: v_pha
%       uv_flag=-1: eta_mean; vn_mean (maxLev*nNodes_obc)
%
%       nf_[mod]=[1,3,4]: nf_temp; nf_salt; nf_[mod]
%       nf_[mod] is of m*n, m means the # of open boundaries, n means the #
%       of tracers in the same module, since the # of tracers is greater
%       than 1 for some modules.
%
%% Notes
% The dimension of 'bc_flags' is M*N, where M depends on the activated
% tracer modules, and M should be at least 4 (purely hydrodynamic case). N
% is the # of open boundaries. N can be 1 in any cases, meaning that the
% treatments of all boundary variables are identical on all boundaries.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 30 Nov 2024.
% Email: wwu@vims.edu
% 
% See also: add_nodal_factors and get_fes2014_tide

%% Parse inputs
if isscalar(find(size(bc_flags)~=1))
    bc_flags = repmat(bc_flags(:)', [Mobj.obc_counts,1]);
end

if size(bc_flags,1)~=1 && size(bc_flags,1)~=Mobj.obc_counts
    error('the # of open boundary segments does not match bc_flags!')
end
if size(bc_flags, 2)<4
    error('the 2nd dimension of bc_flags cannot be less than 4!')
end

if ~isfield(TideForc, 'cutoff_depth')
    error('cutoff_depth is not available!')
end
%% Check the validity of the data
field_list = fieldnames(TideForc);
for ii = 1:numel(field_list)
    tide_var = field_list{ii};
    if strcmpi(tide_var, 'tide_list'); continue; end
    if any(isnan(TideForc.(tide_var)(:)))
        error(['NaNs were found in the variable "', tide_var, '"'])
    end
end
%% Prepare tide info.
tide_list = TideForc.tide_list;
nTides = length(tide_list);

if nTides~=0
    load('tide_fac_constants.mat', 'const');
    tide_pool = cellstr(strtrim(string(const.name)));
    type_pool = const.doodson(:,1); type_pool = fillmissing(type_pool, 'previous'); type_pool = min(type_pool,2);

    tide_types = nan(nTides,1);
    tide_amps = nan(nTides,1);
    tide_freqs = nan(nTides,1);
    for iTide = 1:nTides
        tide_name = strtrim(tide_list{iTide});
        ind_tide = find(strcmpi(tide_pool, tide_name));

        tide_types(iTide) = type_pool(ind_tide,1);
        tide_amps(iTide) = const.potentialamp(ind_tide);
        tide_freqs(iTide) = const.freq(ind_tide)/3600*2*pi;  % rad/s
    end

    TideForc = add_nodal_factors(Mobj, TideForc);

    nodal_factor = TideForc.nodal_factor;
    eq_arg = TideForc.eq_arg;
end

used_mods = check_tracer_module(Mobj);
%% Write tide data
head_line = datestr(now, 'mmm/dd/yyyy HH:MM:SS'); %#ok<TNOW1,DATST>

filepath = fullfile(Mobj.aimpath, 'bctides.in');
fid = fopen(filepath,'w');
fprintf(fid, [head_line, '\n']);

% earth tidal potential and cut-off depth
fprintf(fid, '%d %d.  !# of tidal potential and cut-off depths\n', nTides, TideForc.cutoff_depth);

for iTide = 1:nTides
    % tidal constituent name
    fprintf(fid, [upper(tide_list{iTide}), '\n']);
    % tidal species # (0: declinational; 1: diurnal; 2: semi-diurnal), amplitude constants, angular frequency, nodal factor, earth equilibrium argument (in degrees);
    fprintf(fid, '%d% 8.6f% 12.6e% 8.5f% 9.5f\n', tide_types(iTide), tide_amps(iTide), tide_freqs(iTide), nodal_factor(iTide), eq_arg(iTide));
end

% total # of tidal boundary forcing frequencies
fprintf(fid, [num2str(nTides, '%d'), ' ! total # of tidal boundary forcing frequencies\n']);

for iTide = 1:nTides
    % tidal constituent name
    fprintf(fid, [upper(tide_list{iTide}), '\n']);
    % angular frequency (rad/s), nodal factor, earth equilibrium argument (in degrees) for constituent
    fprintf(fid, '%13.6e% 8.5f% 9.5f\n', tide_freqs(iTide), nodal_factor(iTide), eq_arg(iTide));
end
%% Write B.C for different boundary segments
% # of open boundary segments
fprintf(fid, [num2str(Mobj.obc_counts, '%d'), '  ! # of open boundary segments\n']);

for iSeg = 1: Mobj.obc_counts
    % # of nodes on the open boundary segment j (corresponding to hgrid.gr3), B.C. flags for elevation, velocity, 
    % temperature, salinity, and (optionally) for each tracer module invoked (in the order of GEN, AGE, SED3D, EcoSim, ICM, CoSiNE, FIB, and TIMOR)   
    fprintf(fid, [repmat('%d ', 1, size(bc_flags,2)+1), '! the #',num2str(iSeg), ' open boundary \n'],  Mobj.obc_lens(iSeg), bc_flags(iSeg,:));
   
    elev_flag = bc_flags(iSeg,1);
    uv_flag = bc_flags(iSeg, 2);
    seg_nodes = Mobj.obc_nodes(1:Mobj.obc_lens(iSeg), iSeg);
    ind_nodes = minfind(Mobj.obc_nodes_tot, seg_nodes);
    % ====================== Elevation PART ===========================
    switch elev_flag
        case 0
             disp(['the elevation at boundary (#',num2str(iSeg),') is not specified'])
        case 1 
            disp(['the elevation at boundary (#',num2str(iSeg),') is forced by elev.th (time history)'])
        case 2
            disp(['the elevation at boundary (#',num2str(iSeg),') is forced by a constant elevation'])
            fprintf(fid, '%9.5f\n', TideForc.const_elev(iSeg)); 
        case {3, 5}
            if elev_flag==3
                disp(['the elevation at boundary (#',num2str(iSeg),') is forced by tidal elevation'])
            else
                disp(['the elevation at boundary (#',num2str(iSeg),') is forced by elev2D.th.nc + tidal elevation'])
            end
            for iTide = 1:nTides
                fprintf(fid, [lower(tide_list{iTide}), '\n']);
                ind_tide = find(contains(TideForc.tide_list, tide_list{iTide}));
                elev_amp = TideForc.elev_amp(:,ind_tide);
                elev_pha = TideForc.elev_pha(:,ind_tide);
                for iNodes = ind_nodes(:)'
                    fprintf(fid, '%8.6f% 10.6f\n', elev_amp(iNodes), elev_pha(iNodes));
                end
            end
        case 4
            disp(['the elevation at boundary (#',num2str(iSeg),') is forced by elev2D.th.nc'])
    end
    % ====================== Velocity PART ===========================
    switch uv_flag
        case 0
            if elev_flag==0
                error(['both elevation and velocity at boundary (#',num2str(iSeg),') are not specified!'])
            else
                disp(['the velocity at boundary (#',num2str(iSeg),') is not specified'])
            end
        case 1
            disp(['the velocity at boundary (#',num2str(iSeg),') is forced by flux.th (time history)'])
        case 2
            disp(['the velocity at boundary (#',num2str(iSeg),') is forced by a constant discharge'])
            fprintf(fid, '%9.5f\n', TideForc.const_flow(iSeg));  % constant discharge (note that a negative number means inflow)
        case {3, 5}
            if uv_flag == 3
                disp(['the velocity at boundary (#',num2str(iSeg),') is forced by tidal velocity'])
            else
                disp(['the velocity at boundary (#',num2str(iSeg),') is forced by uv3D.th.nc + tidal velocity'])
            end
            if isfield(TideForc, 'u_amp')
                for iTide = 1:nTides
                    fprintf(fid, [lower(tide_list{iTide}), '\n']);
                    ind_tide = find(contains(TideForc.tide_list, tide_list{iTide}));
                    u_amp = TideForc.u_amp(:,ind_tide);
                    u_pha = TideForc.u_pha(:,ind_tide);
                    v_amp = TideForc.v_amp(:,ind_tide);
                    v_pha = TideForc.v_pha(:,ind_tide);
                    for iNodes = ind_nodes(:)'
                        fprintf(fid, '%8.6f% 10.6f% 8.6f% 10.6f\n', u_amp(iNodes), u_pha(iNodes), v_amp(iNodes), v_pha(iNodes));
                    end
                end
            else
               error('tidal velocity data is not found!')
            end
        case 4
            disp(['the velocity at boundary (#',num2str(iSeg),') is forced by uv3D.th.nc (without nudging)'])
        case -4
            disp(['the velocity at boundary (#',num2str(iSeg),') is forced by uv3D.th.nc (with nudging)'])
            % relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)
            fprintf(fid, '%3.2f %3.2f\n', TideForc.nf_inflow(iSeg), TideForc.nf_outflow(iSeg));
        case -1
            disp('Flather type radiation b.c. (iettype must be 0 in this case)')
            if elev_flag ~= 0
                error('elev_flag must be zero for Flanther radiation case!')
            end
            fprintf(fid, '%5.3f\n', TideForc.eta_mean(iSeg));  % mean elev at each node
            for iNodes = ind_nodes(:)'
                vn_mean = TideForc.vn_mean(:, iNodes);
                vn_mean(isnan(vn_mean)) = []; vn_mean = vn_mean(:)';
                fprintf(fid, [repmat('%5.3f ', 1, size(vn_mean,2)), '\n'], vn_mean);  % mean normal velocity at the node (at all levels)
            end
    end
    % ====================== Tracer PART ===========================
    tracer_mod_list = ['temp', 'salt', used_mods];
    nTracers = numel(tracer_mod_list);
    for ii = 1:nTracers
        tracer_mod_name = tracer_mod_list{ii};
        add_tracer_module(fid, TideForc, Mobj, tracer_mod_name, bc_flags(iSeg, ii+2), iSeg);
    end
end

fclose(fid);
disp('bctides.in has been created successfully!')
end

function used_mods = check_tracer_module(Mobj)

tracer_mod_list = {'gen', 'age', 'sed3d', 'ecosim', 'icm', 'cosine', 'fib', 'timor'};

used_mods = {};
for iMod = 1:8
    tracer_mod_name = tracer_mod_list{iMod};
    if strcmpi(Mobj.(['use_', tracer_mod_name]), 'yes')
        used_mods = [used_mods; tracer_mod_name]; %#ok<AGROW> 
    end
end
end

function add_tracer_module(fid, TideForc, Mobj, tracer_mod_name, mod_flag, iSeg)
% Add tracer boundary condition

tracer_mod_list = {'temp','salt','gen', 'age', 'sed3d', 'ecosim', 'icm', 'cosine', 'fib', 'timor'};
ind_mod = find(contains(tracer_mod_list, tracer_mod_name));
nTracers = Mobj.tracer_counts(ind_mod);

if isfield(TideForc, ['nf_', tracer_mod_name])
    nf_mod = TideForc.(['nf_', tracer_mod_name]);
    if isscalar(nf_mod)
        nf_mod = repmat(nf_mod, [Mobj.obc_counts, nTracers]);
    else
        nf_mod = reshape(nf_mod, [Mobj.obc_counts, nTracers]);
    end
end
if isfield(TideForc, ['const_', tracer_mod_name])
    const_mod = TideForc.(['const_', tracer_mod_name]);
    if isscalar(const_mod)
        const_mod = repmat(const_mod, [Mobj.obc_counts, nTracers]);
    else
        const_mod = reshape(const_mod, [Mobj.obc_counts, nTracers]);
    end
end

switch mod_flag
    case 0
        disp([tracer_mod_name, ' at boundary (#', num2str(iSeg),') is not specified'])
    case 1
        disp([tracer_mod_name, ' at boundary (#', num2str(iSeg),') is forced by *.th'])
        fprintf(fid, [repmat('%3.2f ', 1, size(nf_mod,2)), '\n'], nf_mod(iSeg,:));  % nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow
    case 2
        
        disp([tracer_mod_name, ' at boundary (#', num2str(iSeg),') is forced by a constant value'])
        fprintf(fid, [repmat('%5.3f ', 1, size(const_mod,2)), '\n'], const_mod(iSeg,:));
    case 3
        disp([tracer_mod_name, ' at boundary (#', num2str(iSeg),') is forced by initial profile'])
        fprintf(fid, [repmat('%3.2f ', 1, size(nf_mod,2)), '\n'], nf_mod(iSeg,:));
    case 4
        disp([tracer_mod_name, ' at boundary (#', num2str(iSeg),') is forced by *._3D.th.nc'])
        fprintf(fid, [repmat('%3.2f ', 1, size(nf_mod,2)), '\n'], nf_mod(iSeg,:));
end
end

%% Some Notes
% !...  Count # of tracer models and tracers
% !...  Each trace module has a pre-defined ID as follows:
% !     1: T
% !     2: S
% !     3: GEN
% !     4: AGE
% !     5: SED3D
% !     6: EcoSim
% !     7: ICM
% !     8: CoSINE
% !     9: Feco
% !    10: TIMOR
% !    11: FABM
% !    12: DVD numerical mixing analysis of Klingbeit

% ----------- Elevation B.C. section ------------
% Elev. B.C. = 1, time history of elevation on this boundary
% no input in bctides.in; time history of elevation is read in from elev.th (ASCII);
% Elev. B.C. = 2, this boundary is forced by a constant elevation;
% ethconst !constant elevation value for this segment
% Elev. B.C. = 3, this boundary is forced by tides (amp. and pha.)
% Elev. B.C. = 4, space- and time-varying input
% no input in this file; time history of elevation is read in from elev2D.th.nc (netcdf);
% Elev. B.C. = 5, combination of ‘3’ and ‘4’
% time history of elevation is read in from elev2D.th.nc, and then added to tidal B.C. specified below
% Elev. B.C. = 0
% elevations are not specified for this boundary (in this case the velocity must be specified).

% ----------- Velocity B.C. section ------------
% Velo. B.C. = 1,  time history of discharge on this boundary
% no input in this file; time history of discharge is read in from flux.th (ASCII)
% Velo. B.C. = 2,  this boundary is forced by a constant discharge
% vthconst !constant discharge (note that a negative number means inflow)
% Velo. B.C. = 3,  vel. (not discharge!) is forced in frequency domain
% tidal amplitude and phase for (u,v) at each node on this open boundary
% Velo. B.C. = 4, time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)
% Velo. B.C. = -4, time history of velocity (not discharge!) is read in
% from uv3D.th.nc (netcdf). In addition, relaxation constants for inflow and
% outflow (between 0 and 1 with 1 being strongest nudging
% Velo. B.C. = 5, combination of '4' and '3'
% time history of velocity (not discharge!) is read in from uv3D.th.nc
% (netcdf) and then added to tidal velocity specified below
% Velo. B.C. = -1, Flanther type radiation b.c. (iettype must be 0 in this case)
% Velo. B.C. = 0, vel. not specified, no input needed

% ----------- Temperature B.C. section ------------
% Temp. B.C. = 1, time history of temperature on this boundary
% tobc ! nudging factor (between 0 and 1 with 1 being strongest nudging) for
% inflow; time history of temperature will be read in from TEM_1.th (ASCII)
% Temp. B.C. = 0, temperature not specified
% Temp. B.C. = 2, this boundary is forced by a constant temperature
% tthconst !constant temperature on this segment
% tobc !nudging factor (between 0 and 1) for inflow
% Temp. B.C. = 3, initial temperature profile for inflow
% tobc !nudging factor (between 0 and 1) for inflow
% Temp. B.C. = 4, 3D input
% tobc !nudging factor (between 0 and 1); time history of temperature is
% read in from TEM_3D.th.nc (netcdf).

% ----------- Salinity B.C. section ------------
% Salinity B.C. section is similiar to the Temperature Part.

% The cut-off depth is used to save a little computational time as tidal potential
% is negligible in shallow water. To include all depths in calculation, set
% a negative depth like -100.

% ELEV, VEL, TEM, SAL, GEN, AGE, SED3D, EcoSim, ICM, CoSiNE, FIB, and TIMOR 

% tide_pool = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
% tide_types = [2 2 2 2 1 1 1 1];
% tide_amps = [0.112841 0.242334 0.046398 0.030704 0.141565 0.046843 0.100514 0.019256];
% tide_freqs = [0.145444E-03 0.140519E-03 0.137880E-03 0.145842E-03 0.729212E-04 0.725229E-04 0.675977E-04 0.649585E-04];








