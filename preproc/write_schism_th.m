function write_schism_th(Mobj, suffixName, varMatrix, time_steps)
% 
% 
%% Syntax
% 
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
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2022-06-06.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
% This includes elev.th, flux.th, TEM_1.th, SAL_1.th etc
% varMatrix (nTimes*nNodes)
% the time interval of the provided 'timeSteps' must be larger than Mobj.dt

%%
if isnan(sum(varMatrix(:)))
    error('NaN values were mixed into varMatrix, please check!')
end

if size(varMatrix,1) ~= length(time_steps)
    varMatrix = varMatrix';
end

varNew = multi_interp1(time_steps, varMatrix, Mobj.time, 1);
timeTicks = seconds(Mobj.time - Mobj.time(1));

fileName = [Mobj.aimpath, suffixName, '.th'];
fid = fopen(fileName,'wt');
formatStr = ['%d.', repmat('% 11.3f', 1, size(varNew,2)), '\n'];
for iTime = 1:length(timeTicks)
    progressbar(iTime/length(timeTicks))
    fprintf(fid, formatStr, timeTicks(iTime), varNew(iTime,:));
end
fclose(fid);
disp([suffixName, '.th has been created successfully!'])
end