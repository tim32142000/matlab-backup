clear all
[t,yy]=ode45(@damped1,[0 30],[0.5 0.8]); 
figure,plot(t,yy(:,1)), hold on, plot(t,yy(:,2),'r')