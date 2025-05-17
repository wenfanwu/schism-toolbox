function config = read_schism_nml(filepath)
% Read the *.nml file for SCHISM.
%
%% Syntax
% config = read_schism_nml(filepath)
%
%% Description
% config = read_schism_nml(filepath) reads the parameters in *.nml file.
%
%% Examples 
% filepath = 'schism-master\sample_inputs\icm.nml';
% config = read_schism_nml(filepath);
%
%% Input Arguments
% filepath - file path; char
%       the absolute file path of *nml file.
%
%% Output Arguments
% config - model config; datastruct
%       the datastruct used to store all parameters.
%
%% Notes
% This function was created with the help of ChatGPT.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 16 May 2025. 
% Email: wwu@vims.edu
% 
% See also: read_schism_gr3

%% Parse inputs
fid = fopen(filepath, 'r');
if fid == -1
    error('Cannot open file: %s', filepath);
end

config = struct();
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '!') || startsWith(line, '&') || strcmp(line, '/')
        continue;
    end
    line = strtok(line, '!');

    if contains(line, '=')
        [raw_key, val] = strtok(line, '=');
        val = strtrim(val(2:end));  % value string
        raw_key = strtrim(raw_key);

        % Check if the key is indexed (e.g., var(1))
        idx_match = regexp(raw_key, '(\w+)\((\d+)\)', 'tokens');
        if ~isempty(idx_match)
            key = idx_match{1}{1};
            idx = str2double(idx_match{1}{2});
        else
            key = raw_key;
            idx = [];
        end

        % Parse value
        if startsWith(val, '''') && endsWith(val, '''')
            parsed_value = strrep(val, '''', '');
        else
            % handle Fortran style：'1.d0' → '1.0d0'，d/D →  e
            val = regexprep(val, '(?<=\d)\.(?=d[\+\-\d])', '.0');
            val = regexprep(val, '[dD]', 'e');

            parts = regexp(val, '[^,\s]+', 'match');
            nums = str2double(parts);
            if all(~isnan(nums))
                if isscalar(nums)
                    parsed_value = nums(1);
                else
                    parsed_value = nums(:)';
                end
            else
                parsed_value = val;
            end
        end

        % Store to config
        if isempty(idx)
            config.(key) = parsed_value;
        else
            if ~isfield(config, key)
                config.(key) = [];
            end
            config.(key)(idx) = parsed_value;
        end
    end
end
fclose(fid);
end
