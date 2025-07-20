% 
%  T+-epsilon control on Chaotic Compass
%  time is in units of 1/(the driving frequency omega)
%  d^2 theta/dt^2+ gam* d theta/dt +b1 cos(theta)+b2 sin(theta)*cos(2\pi*t)=0
% ******************************************************************************
% PRO compass
clear all
nwindow=256;
dt = 0.001;
t = 0.0;
totaltime=1600;
nstep=(totaltime/dt);
istep=1/dt; % ~ to period of external force
iton=200/dt;
itoff=2.*nstep;
   i3=1750/dt;
eps0=0.5;
% epstilde=eps0*0.5;
eps=0.0;
b1=36.  ; % for cos(twopi*t) driving:  b1=36, b2=98, gam=6 gives period-2.  b1=36, b2=99.5 gam=6 gives period-4
b20=99.;    %                     b2=97.5-- period-1   b2=97.6---period-2  b2=101 -- chaotic  b2=105 --period-3
gam=6.;
alpha=1.; % w0 in the equation, frequency of external force

% ;gam=2.2;
% ;b1=0.32;
% ;b2=11.24;



dtg=dt*gam;
dt2=dt*dt;
twopi=6.283185307;
strength=0.;
% main loop starts from here
%x(1)=0.5;
x(1)=30.0/180*pi;
omega(1)=0.;% omega is d theta/dt
x(2)=x(1)+omega(1)*dt;
omega(2)=omega(1);
time(1)=0.;
time(2)=dt;
jj=1;
b2=b20;
xxn(1)=x(1);
for it=3:nstep
t=it*dt;
 if(it > iton && it < itoff), eps=2.; end
%     ;        if(it ge 350 and it le itoff) then eps=eps0/3
if(jj > 3) 
    if(xxn(jj-1) < xxn(jj-2)), b2=b20+eps; end
    if(xxn(jj-1) > xxn(jj-2)), b2=b20-eps; end
end
 x(it)=((2+dtg)*x(it-1)-x(it-2)-dt2*(b1*cos(x(it-1))+1*b2*sin(x(it-1))*cos(1*alpha*twopi*t)))/(1+dtg);
%   x[it]=((2-dtg)*x[it-1]+(dtg-1)*x[it-2]-dt2*(b1*cos(x[it-1])+b2*sin(x[it-1])*cos(alpha*twopi*t)))/(1+0)
omega(it)=(x(it)-x(it-1))/dt;
time(it)=t;
 if( mod(it,istep)== 0)
 xxn(jj)=x(it);
 jj=jj+1;
 end
end % it loop

figure,plot(time(9*nstep/10:nstep),x(9*nstep/10:nstep)),xlabel('t'), ylabel('theta(t)')
title('time series')
figure,plot(x(nstep/2:nstep),omega(nstep/2:nstep)),xlabel ('\theta'), ylabel('d\theta/dt') 
title('phase space')

for j=1:totaltime 
it=j/dt;
xn(j)=x(it-00);
omegan(j)=omega(it-00);
end
figure,plot( xxn(iton*dt+5-1:totaltime-1),xxn(iton*dt+5:totaltime),'.','MarkerSize',6),xlabel ('\theta n'), ylabel('\theta_{n+1}') 
title('shape of chaotic attractor')
figure,plot(xxn,'.-')%,ylim([-3.,-2.70]),
xlabel('n'), ylabel('\theta n') 
title('state of system like map model')
%figure,plot(xxn(25:150), xxn(26:151)),xlabel('\theta_n'), ylabel('\theta n+1') 

