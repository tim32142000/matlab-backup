function dy = damped1(t,y)
w = 2*pi ;
gamma = 0.01 ;
dy = zeros(2,1) ;
dy(1) = y(2) ;
dy(2) = -gamma*y(2) - w*w*y(1) ;