clear all
[t,yy]=ode45(@harmonic1,[0 10],[1. 1.5]); % tspan=[0 10]; intial=[1. 1.5]
figure,plot(t,yy(:,1)), hold on, plot(t,yy(:,2),'r')