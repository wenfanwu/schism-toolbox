function write_schism_bctides(Mobj, TideForc, bc_flags)
% Write the vgrid.in file for SCHISM.
% 
%% Syntax
% write_schism_bctides(Mobj, TideForc, bc_flags)
%
%% Description 
% write_schism_bctides(Mobj, TideForc, bc_flags) creates bctide.in file
%
%% Examples
% tideList = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
% TideForc = get_fes2014_tide(Mobj, tideList);   
% write_schism_bctides(Mobj, TideForc, [5 5 4 4])
%
%% Input Arguments
% Mobj --- the mesh object
% TideForc --- datastruct containing the tidal forcing data, which can be
% obtained from the 'get_fes2014_tide.m' function.
% bc_flags --- open boundary types specified for different tracers. 
%
%% Output Arguments
% None
% 
%% Notes
% The dimension of 'bc_flags' is M*N, where M depends on the activated
% tracer modules, and M should be at least 4 (purely hydrodynamic case). N
% is the # of open boundaries. N can be 1 in any cases, meaning that the
% treatments of all boundary variables are identical on all boundaries.
% Instead, you can define different flags for different boundaries.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: add_nodal_factors and get_fes2014_tide

%% Parse inputs
if size(bc_flags,1)~=1 && size(bc_flags,1)~=Mobj.obc_counts
    error('The # of B.C. segments does not correspond to the second dimension of bc_flags')
end

used_mods = check_tracer_module(Mobj);

ntrs = numel(used_mods);  % not strict, as the tracers can exceed the number of module
if ntrs~=0
    err_mods = '';
    for ii = 1:ntrs
        err_mods = [err_mods, '+', used_mods{ii}]; %#ok<AGROW> 
    end
    error(['please specify B.C. flags for the tracer modules (',err_mods(2:end),')'])
end

%% Define some default values
if isfield(TideForc, 'cutoff_depth')
    cutoff_depth = TideForc.cutoff_depth;
else
    cutoff_depth = -20;
end

tideList = TideForc.tide_list;
tideNums = length(tideList);
tideNumsPot = tideNums;

% this part needs improvement in the future
tideTypes = [2 2 2 2 1 1 1 1];
tideAmps = [0.112841 0.242334 0.046398 0.030704 0.141565 0.046843 0.100514 0.019256];
tideFreqs = [0.145444E-03 0.140519E-03 0.137880E-03 0.145842E-03 0.729212E-04 0.725229E-04 0.675977E-04 0.649585E-04];

nodalFactors = TideForc.nf;
eqArg = TideForc.eq_arg;

headLine = datestr(Mobj.time(1), 'mm/dd/yyyy HH:MM:SS UTC');

fileName = [Mobj.aimpath, 'bctides.in'];
fid = fopen(fileName,'wt');
fprintf(fid, [headLine, '\n']);

%% Begin to write
% Earth tidal potential and Cut-off depth
fprintf(fid, '%d %d.  !# of tidal potential and cut-off depths\n', tideNumsPot, cutoff_depth);
if tideNumsPot ~= 0
    for iTide = 1:tideNumsPot
        fprintf(fid, [tideList{iTide}, '\n']);
        fprintf(fid, '%d% 8.6f% 12.6e% 8.5f% 9.5f\n', tideTypes(iTide), tideAmps(iTide), tideFreqs(iTide),nodalFactors(iTide), eqArg(iTide));
    end
end
% Activated tidal consitituents
fprintf(fid, [num2str(tideNums, '%d'), ' ! total # of tidal boundary forcing frequencies\n']);
for iTide = 1:tideNums
    fprintf(fid, [tideList{iTide}, '\n']);
    fprintf(fid, '%13.6e% 8.5f% 9.5f\n', tideFreqs(iTide),nodalFactors(iTide), eqArg(iTide));
end

%% B.C for different boundary segments
if size(bc_flags,1) ~= Mobj.obc_counts
    bc_flags = repmat(bc_flags(:)', Mobj.obc_counts, 1);
end
fprintf(fid, [num2str(Mobj.obc_counts, '%d'), '  ! # of open boundary segments\n']);

for iSeg = 1: Mobj.obc_counts
    fprintf(fid, [repmat('%d ', 1, size(bc_flags,2)+1), '\n'],  Mobj.obc_lens(iSeg), bc_flags(iSeg,:));
   
    elev_flag = bc_flags(iSeg,1);
    uv_flag = bc_flags(iSeg, 2);
    segNodes = Mobj.obc_nodes(1:Mobj.obc_lens(iSeg),iSeg);
    indNodes = minfind(Mobj.obc_nodes_tot, segNodes);
    % ====================== Elevation PART ===========================
    switch elev_flag
        case 1
            disp('time history of elevation is read in from elev.th (ASCII)')
        case 2
            fprintf(fid, '%9.5f\n', TideForc.const_elev);
        case {3, 5}
            for iTide = 1:tideNums
                fprintf(fid, [lower(tideList{iTide}), '\n']);
                indTide = find(contains(TideForc.tide_list, tideList{iTide}));
                elev_amp = TideForc.elev_amp(:,indTide);
                elev_pha = TideForc.elev_pha(:,indTide);
                for iNodes = indNodes(:)'
                    fprintf(fid, '%8.6f% 10.6f\n', elev_amp(iNodes), elev_pha(iNodes));
                end
            end
        case 4
            disp('time history of elevation is read in from elev2D.th.nc')
    end
    % ====================== Velocity PART ===========================
    switch uv_flag
        case 0
            disp('B.C. velocity is not specified! No input is needed in the bctides.in')
        case 1
            disp('time history of discharge is given on this boundary')
        case 2
            disp('this boundary is forced by a constant discharge')
            fprintf(fid, '%9.5f\n', TideForc.const_flow);
        case {3, 5}
            if uv_flag == 3
                disp('vel. (not discharge!) is forced in frequency domain')
            end
            if uv_flag == 5
                disp('time history of velocity (not discharge!) is read in from uv3D.th.nc and then added to the tidal velocity')
            end
            if isfield(TideForc, 'u_amp')
                for iTide = 1:tideNums
                    fprintf(fid, [lower(tideList{iTide}), '\n']);
                    indTide = find(contains(TideForc.tide_list, tideList{iTide}));
                    u_amp = TideForc.u_amp(:,indTide);
                    u_pha = TideForc.u_pha(:,indTide);
                    v_amp = TideForc.v_amp(:,indTide);
                    v_pha = TideForc.v_pha(:,indTide);
                    for iNodes = indNodes(:)'
                        fprintf(fid, '%8.6f% 10.6f% 8.6f% 10.6f\n', u_amp(iNodes), u_pha(iNodes), v_amp(iNodes), v_pha(iNodes));
                    end
                end
            else
                warning('velocity part is not found in the given struct!')
            end
        case 4
            disp('time history of velocity (not discharge!) is read in from uv3D.th.nc')
        case -4
            disp('time history of velocity (not discharge!) is read in from uv3D.th.nc')
            %relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)
            fprintf(fid, '%5.3f %5.3f\n', TideForc.relax_const_in, TideForc.relax_const_out);
        case -1
            disp('Flanther type radiation b.c. (iettype must be 0 in this case)')
            if elev_flag ~= 0
                error('The elev_flag must be zero for Flanther radiation case!')
            end
            fprintf(fid, '%d\n', TideForc.eta_mean);
    end
    % ====================== Tracer PART ===========================
    tracer_list = ['temp', 'salt', used_mods];
    nTracers = numel(tracer_list);
    for ii = 1:nTracers
        tracer_mod_name = tracer_list{ii};
        add_tracer_module(fid, TideForc, tracer_mod_name, bc_flags(iSeg, ii+2));
    end
end

% --------------- Add River Flux at the end, if necessary
% if river_flag == 1
%     for iRiver = 1:RivFlux.riverNums
%         riverName = RivFlux.riverNames{RivFlux.ind(iRiver)};
%         fprintf(fid, [repmat('%d ', 1, size(RivFlux.bcFlags,2)+1), '  !',riverName, '\n'],  RivFlux.branchLens(iRiver), RivFlux.bcFlags(iRiver,:));
%         
%         elev_flag = RivFlux.bcFlags(iRiver,1);
%         uv_flag = RivFlux.bcFlags(iRiver,2);
%         temp_flag = RivFlux.bcFlags(iRiver,3);
%         salt_flag = RivFlux.bcFlags(iRiver,4);
%         
%         indNodes = RivFlux.riverNodes(iRiver:iRiver+RivFlux.branchLens(iRiver)-1);
%         
%         switch uv_flag
%             case 1
%                 disp('time history of discharge is given on this boundary')
%         end
%         switch temp_flag
%             case 1
%                 disp('time history of temperature on this boundary; time history of temperature will be read in from TEM_1.th')
%                 fprintf(fid, '%2.1f\n', RivFlux.nf_temp);    %nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow
%         end
%         switch salt_flag
%             case 1
%                 disp('time history of saltlinity on this boundary; time history of salinity will be read in from SAL_1.th')
%                 fprintf(fid, '%2.1f\n', RivFlux.nf_salt);    %nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow
%         end
%     end
% end

fclose(fid);

disp('bctides.in has been created successfully!')
end

function used_mods = check_tracer_module(Mobj)
mod_list = {'gen', 'age', 'sed3d', 'ecosim', 'icm', 'cosine', 'fib', 'timor'};
used_mods = {};
for iMod = 1:8
    mod_name = mod_list{iMod};
    if strcmpi(Mobj.(['use_', mod_name]), 'yes')
        used_mods = [used_mods; mod_name]; %#ok<AGROW> 
    end
end
end

function add_tracer_module(fid, TideForc, tracer_mod_name, mod_flag)

abbr_name = upper(tracer_mod_name(1:3));
switch mod_flag
    case 0
        disp([tracer_mod_name, ' is not specified'])
    case 1
        nf_mod = TideForc.(['nf_', tracer_mod_name]);
        disp(['time history of salinity on this boundary; time history of temperature will be read in from ', abbr_name,'.th'])
        fprintf(fid, '%d\n', nf_mod);    % nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow
    case 2
        const_mod = TideForc.(['const_', tracer_mod_name]);
        disp(['this boundary is forced by a constant value for ', abbr_name])
        fprintf(fid, '%5.3f\n', const_mod);
    case 3
        nf_mod = TideForc.(['nf_', tracer_mod_name]);
        disp('initial salinity profile for inflow')
        fprintf(fid, '%d\n', nf_mod);
    case 4
        nf_mod = TideForc.(['nf_', tracer_mod_name]);
        disp(['time history of salinity is read in from ', abbr_name,'_3D.th.nc'])
        fprintf(fid, '%d\n', nf_mod);
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








