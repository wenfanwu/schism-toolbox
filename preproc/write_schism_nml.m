function write_schism_nml(filepath, config)
% Write*.nml file for SCHISM (no comments).
%
%% Syntax
% write_schism_nml(filepath, config)
%
%% Description
% write_schism_nml(filepath, config) write nml file for SCHISM.
%
%% Examples 
% config = read_schism_nml('schism-master\sample_inputs\param.nml');
% config.dt = 120;
% write_schism_nml('schism-master\sample_inputs\param2.nml', config);
%
%% Input Arguments
% filepath - file path; char
%       the absolute file path of *nml file.
% config - model config; datastruct
%       the datastruct used to store all parameters. config is usually
%       obtained from 'read_schism_nml'.
%
%% Output Arguments
% None
%
%% Notes
% This function was created with the help of ChatGPT
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 16 May 2025. 
% Email: wwu@vims.edu
% 
% See also: read_schism_nml

%% Parse inputs
% Variables to be written in expanded form with indexing
expand_vars = {'iof_hydro', 'iof_wwm', 'iof_gen', 'iof_sed', 'iof_sed2d', 'iof_ice', 'iof_ana', ...
    'iadjust_mass_consv0', 'inu_tr', 'lev_tr_source', 'flag_ic', 'stemp_dz'};

fid = fopen(filepath, 'w');
if fid == -1
    error('Cannot open file for writing: %s', filepath);
end

% Start of namelist block
fprintf(fid, '&SCHISM_NAMELIST\n');

vars = fieldnames(config);
for i = 1:length(vars)
    key = vars{i};
    val = config.(key);

    % Case 1: Expand indexed format (one line per element)
    if ismember(key, expand_vars) && isnumeric(val) && isvector(val)
        for j = 1:length(val)
            if mod(val(j), 1) == 0
                fprintf(fid, '  %s(%d) = %d\n', key, j, val(j));
            else
                fprintf(fid, '  %s(%d) = %.10g\n', key, j, val(j));
            end
        end

        % Case 2: Scalar or array in single line with 3-space separators
    else
        if ischar(val)
            fprintf(fid, '  %s = ''%s''\n', key, val);
        elseif isnumeric(val)
            if isscalar(val)
                if mod(val, 1) == 0
                    fprintf(fid, '  %s = %d\n', key, val);
                else
                    fprintf(fid, '  %s = %.10g\n', key, val);
                end
            else
                if all(mod(val, 1) == 0)
                    str_vals = arrayfun(@(x) sprintf('%d', x), val, 'UniformOutput', false);
                else
                    str_vals = arrayfun(@(x) sprintf('%.10g', x), val, 'UniformOutput', false);
                end
                joined = strjoin(str_vals, '   ');
                fprintf(fid, '  %s = %s\n', key, joined);
            end
        else
            warning('Skipping unsupported field: %s', key);
        end
    end
end

% End of namelist block
fprintf(fid, '/\n');
fclose(fid);
end



