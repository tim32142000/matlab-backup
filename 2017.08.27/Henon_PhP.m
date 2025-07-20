% Matlab code (see [4]) to demonstrate H?non map strange attractor
clc; clear all;
% define the parameters
%a=input('a = ');
a=1.3;
%b=input('b = ');
b=0.3;
% specify the initial conditions.
%x0=input('x0 = ');
x0=rand(1)*3-1.5 ;
%y0=input('y0 = ');
y0=rand(1)*3-1.5 ;
%n=input('Maximum number of iterations = ');
n=2000;
x=zeros(1,n+1);
y=zeros(1,n+1);
x(1)=x0;
y(1)=y0;
% main routine
for i=1:n         % iterates the H?nonmap
    x(i+1)=1-a*(x(i)^2)+y(i);
    y(i+1)=b*x(i);
end
plot(x(1001:2000),y(1001:2000),'.k','LineWidth',.5,'MarkerSize',5);
xlabel('x_n'); ylabel('y_n');
title(['Henon map: a= ',num2str(a),', b= ',num2str(b),', (x_0,y_0)=(',num2str(x0),',',num2str(y0),')']);
grid
zoom