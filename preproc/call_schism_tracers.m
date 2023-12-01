function Mobj = call_schism_tracers(Mobj)
% Load the Info. of activated tracers
% 
%% Syntax
% Mobj = call_schism_tracers(Mobj)
%
%% Description 
%  Mobj = call_schism_tracers(Mobj) adds the info. of activated modules
%  into Mobj.
%
%% Example
% Mobj = call_schism_tracers(Mobj)
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
% Last Updated on 2022-10-10.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 
%% All modules are off by default
module_list = {'use_gen', 'use_age', 'use_sed3d', 'use_ecosim', 'use_icm', 'use_icm_ph', 'use_cosine', 'use_fib', 'use_timor', ...
    'use_wwm', 'use_ice'};  % non-tracer modules
for imod = 1:numel(module_list)
    mod_name = module_list{imod};
    if ~isfield(Mobj, mod_name)
        Mobj.(mod_name) = 'no';
    else
        Mobj = rmfield(Mobj, mod_name);
        Mobj.(mod_name) = 'yes';
    end
end
%% TEM&SAL (index=1&2)
active_mods = cell(10,1);
active_tracers = {'temp', 'salt'};
tracer_counts = zeros(1,10);

nRows = 100; nMods = 10;
tracer_sheet = cell(nRows, nMods);

tracer_counts(1:2) = 1;
active_mods{1} = 'TEM';
active_mods{2} = 'SAL';

tracer_sheet(1, 1) = {'TEM'};
tracer_sheet(1, 2) = {'SAL'};
tracer_sheet(1, :) = {'TEM', 'SAL', 'GEN', 'AGE', 'SED', 'ECO', 'ICM', 'COS', 'FIB', 'TMR'};

tracer_sheet(2,:) = {'off'};
tracer_sheet(2, 1:2) = {'ON'};


%% USE_GEN (index=3)
if strcmpi(Mobj.use_gen, 'yes')
    active_mods{3} = 'GEN';
    tracer_sheet(2, 3) = {'ON'};

    if isfield(Mobj, 'ntracer_gen')
        if Mobj.ntracer_gen<0
            error('the attibute "ntracer_gen" must be positive integer!')
        end
        tracer_counts(3) = Mobj.ntracer_gen;
    else
        error('the attibute "ntracer_gen" is not specified!')
    end
    if tracer_counts(3)<=0
        error('INIT: ntrs(3)<=0')
    end
    used_tracers = arrayfun(@(x) ['gen_', num2str(x)], 1:Mobj.ntracer_gen, 'UniformOutput', false);
    active_tracers = [active_tracers(:); used_tracers(:)];
    tracer_sheet(3:Mobj.ntracer_gen+2, 3) = used_tracers;
end
%% USE_AGE (index=4)            
if strcmpi(Mobj.use_age, 'yes')
    active_mods{4} = 'AGE'; 
    tracer_sheet(2, 4) = {'ON'};

    if isfield(Mobj, 'ntracer_age')
        if Mobj.ntracer_age<0
            error('the attibute "ntracer_age" must be positive integer!')
        end
        tracer_counts(4) = Mobj.ntracer_age;
    else
        error('please specify the ntracer_age')
    end
    if tracer_counts(4)<=0 || mod(tracer_counts(4),2)~=0
        error('INIT: age requires even # of tracers')
    end
end
%% USE_SED (index=5)
if strcmpi(Mobj.use_sed3d, 'yes')
    active_mods{5} = 'SED3D';
    tracer_sheet(2, 5) = {'ON'};

    if isfield(Mobj, 'sed_class')
        if Mobj.sed_class<0
            error('the attibute "sed_class" must be positive integer!')
        end
        tracer_counts(5) = Mobj.sed_class;  %
    else
        error('the attibute "sed_class" is not specified!')
    end
    if(tracer_counts(5)<=0)
        error('INIT: ntrs(5)<=0')
    end
    used_tracers = arrayfun(@(x) ['sed_', num2str(x)], 1:Mobj.sed_class, 'UniformOutput', false);
    active_tracers = [active_tracers(:); used_tracers(:)];
    tracer_sheet(3:Mobj.sed_class+2, 5) = used_tracers;
end
%% USE_ECO (index=6) 
if strcmpi(Mobj.use_ecosim, 'yes')
    active_mods{6} = 'ECO'; 
    tracer_sheet(2, 6) = {'ON'};

    if isfield(Mobj, 'eco_class')
        if Mobj.eco_class<0
            error('the attibute "eco_class" must be positive integer!')
        end
        tracer_counts(6) = Mobj.eco_class;
    else
        error('please specify the eco_class')
    end
    if(tracer_counts(6)<=0)
        error('INIT: ntrs(6)<=0')
    end
end
%% USE_ICM (index=7)
if strcmpi(Mobj.use_icm, 'yes') || strcmpi(Mobj.use_icm_ph, 'yes')
    active_mods{7} = 'ICM';
    tracer_sheet(2, 7) = {'ON'};

    tracer_counts(7) = 21;
    used_tracers = {'ZB1'; 'ZB2'; 'PB1'; 'PB2'; 'PB3'; 'RPOC'; 'LPOC'; 'DOC'; 'RPON'; 'LPON'; 'DON'; 'NH4'; 'NO3'; 'RPOP'; 'LPOP'; 'DOP'; 'PO4'; 'SU'; 'SA'; 'COD'; 'DO'};
    tracer_sheet(3:23,7) = used_tracers;
    if strcmpi(Mobj.use_icm_ph, 'yes')
        active_mods{7} = 'ICM_pH';
        tracer_sheet{1, 7} = 'ICM_pH';
        
        tracer_counts(7) = 25;
        used_tracers = {'ZB1'; 'ZB2'; 'PB1'; 'PB2'; 'PB3'; 'RPOC'; 'LPOC'; 'DOC'; 'RPON'; 'LPON'; 'DON'; 'NH4'; 'NO3'; 'RPOP'; 'LPOP'; 'DOP'; 'PO4'; 'SU'; 'SA'; 'COD'; 'DO'; ...
            'TIC'; 'ALK'; 'CA'; 'CACO3'}; % additional variables of pH module 
        tracer_sheet(3:27,7) = used_tracers;
    end
    active_tracers = [active_tracers(:); used_tracers(:)];
end
%% USE_COS (index=8)
if strcmpi(Mobj.use_cosine, 'yes')
    active_mods{8} = 'COS';
    tracer_sheet(2, 8) = {'ON'};

    tracer_counts(8) = 13;
    used_tracers = {'NO3'; 'SiO4'; 'NH4'; 'S1'; 'S2'; 'Z1'; 'Z2'; 'DN'; 'DSi'; 'PO4'; 'DOX'; 'CO2'; 'ALK'}';
    active_tracers = [active_tracers(:); used_tracers(:)];
    tracer_sheet(3:15, 8) = used_tracers;
end
%% USE_FIB  (index=9)
if strcmpi(Mobj.use_fib, 'yes')
    active_mods{9} = 'FIB'; % not completed yet
    tracer_sheet(2, 9) = {'ON'};

    tracer_counts(9)=2;
end
%% USE_TIMOR  (index=10)
if strcmpi(Mobj.use_timor, 'yes')
    active_mods{10} = 'TMR'; % not completed yet
    tracer_sheet(2, 10) = {'ON'};
end
%% Tracer Info
tracer_sheet = tracer_sheet(~all(cellfun('isempty', tracer_sheet), 2), :);
tracer_sheet(cellfun('isempty', tracer_sheet)) = {'---'};
active_mods(cellfun('isempty', active_mods)) = [];

Mobj.tracer_sheet = tracer_sheet;
Mobj.active_mods = active_mods;
Mobj.active_tracers = active_tracers;
Mobj.tracer_counts = tracer_counts;
Mobj.nTracers = sum(Mobj.tracer_counts);

disp(['Activated Modules: ', strjoin(active_mods, ', ')])
disp(['A total of ', num2str(sum(tracer_counts)), ' tracers are considered'])
end
























