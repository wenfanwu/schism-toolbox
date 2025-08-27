function sigma = vector_angle(A, B)
% Calculate the anti-clockwise angle from vector A to B.

x1 = A(:,1); y1 = A(:,2);
x2 = B(:,1); y2 = B(:,2);

dot_val = x1.*x2+y1.*y2;
cross_val = x1.*y2 - y1.*x2;
sigma = atan2(cross_val, dot_val)*180/pi;

end