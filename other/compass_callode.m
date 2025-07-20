B1 = 36 ;
B2 = 105 ;

[t,yy]=ode45(@compass_odefn,[0 2001],[0 0.0],[],B1,B2);
%plot(yy(7000:,1),yy(7000:,2) );