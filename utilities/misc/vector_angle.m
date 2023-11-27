function sigma = vector_angle(A, B)
% Calculate the anti-clockwise angle from vector A to B.
%
%% Syntax
% sigma = vector_angle(A, B)
%
%% Description
% sigma = vector_angle(A, B) returns the anti-clockwise angle from vector A
% to B in degrees.
%
%% Example
% A = [1 0];
% B = [1 sqrt(3)];
% sigma = vector_angle(A, B)
% sigma = vector_angle(B, A)
% 
%% Input Arguments
% A --- a matrix of N vectors (N*2), N = 1, 2, 3 ...
% B --- a matrix of N vectors (N*2), N = 1, 2, 3 ...; the size of B and A must be
% consistent.
% 
%% Input Arguments
% sigma --- a vector in N columns, containing the anti-clockwise angle from
% vector A to B.
% 
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 21 Feb. 2022. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: atan2

%% Parse inputs
x1 = A(:,1);
y1 = A(:,2);
x2 = B(:,1);
y2 = B(:,2);

dot_val = x1.*x2+y1.*y2;
cross_val = x1.*y2 - y1.*x2;
sigma = atan2(cross_val, dot_val)*180/pi;

end