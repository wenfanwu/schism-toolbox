function write_schism_source(Mobj, D, tracer_list)
% Write source_sink.in and relevant files for SCHISM
%
%% Syntax
% write_schism_source(Mobj, D, tracer_list)
%
%% Description 
% 
%
%% Examples
%
%
%% Input Arguments
%
%
%% Output Arguments
% 
% 
%% Notes
% the # of sink/source can not be zeros in SCHISM, but you can set the
% flow volumn to be zero to if the sink/source point is not needed. In
% addition, the source/sink values must be single-precision.
% Source/sinks can be specified at an elem more than once, and the code
% will accumulate the volumes, hence if you assign multiple elems for one
% river, make sure that the sum of the volumes is equal to the runoff.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 27 Nov 2024
% Email: wwu@vims.edu
% 
% See also: write_schism_source_nc

%% Parse inputs
if nargin < 3
    tracer_list = {'temp', 'salt'}; 
end

n_src = numel(D.source_elems);
n_snk = numel(D.sink_elems);
if n_src == 0
    D.vsource.flow = zeros(size(D.vsource.time));
    D.msource.temp = zeros(size(D.vsource.time));
    D.msource.salt = zeros(size(D.vsource.time));
end
if n_snk == 0
    D.vsink.flow = zeros(size(D.vsink.time));
    D.sink_elems = 1;
end

%% source_sink.in
fileName = [Mobj.aimpath, 'source_sink.in']; %source_sink.in; specify the source/sink elements
fid = fopen(fileName,'wt');
fprintf(fid, [num2str(max(1,n_src)), '\n']);
for ii = 1:max(1,n_src)
    fprintf(fid, [num2str(D.source_elems(ii)), '\n']);
end
fprintf(fid, '         \n');  % this blank line is required!
fprintf(fid, [num2str(max(1,n_snk)), '\n']);
for ii = 1:max(1,n_snk)
    fprintf(fid, [num2str(D.sink_elems(ii)), '\n']);
end
fclose(fid);

%% vsource.th & vsink.th
write_schism_th(Mobj, 'vsource', abs(D.vsource.flow), D.vsource.time)
write_schism_th(Mobj, 'vsink', abs(D.vsink.flow), D.vsink.time)

%% msource.th
nTracers = numel(tracer_list);
nTimes = length(D.msource.time);
varAll = [];
for ii = 1:nTracers
    tracer_name = tracer_list{ii};
    varTmp = D.msource.(tracer_name);
    if size(varTmp,1) ~= nTimes
        varTmp = varTmp';
    end
    varAll = [varAll varTmp]; %#ok<*AGROW>
end
varAll(isnan(varAll)) = -9999;   % junk values, will not be used

write_schism_th(Mobj, 'msource', varAll, D.msource.time)
end
















