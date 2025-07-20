function dy=harmonic1(t,y);
w=1*2*pi; % f=w/(2pi)
dy = zeros(2,1);    % a column vector
dy(1) = y(2);
dy(2) = -w^2*y(1);