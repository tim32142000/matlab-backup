function dy = compass_odefn(t,y,B1,B2)
w = 2*pi ;
gamma = 6 ;
dy = zeros(2,1) ;
dy(1) = y(2) ;
dy(2) = -gamma*y(2) - B1*sin(y(1) ) + B2*cos(y(1) )*cos(w*t) ;