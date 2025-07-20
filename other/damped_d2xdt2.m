clear all
tau = 0.0001 ;
tot_time = 30 ;
tot_ts = tot_time/tau ;
w = 2*pi ;
gamma = 0.01 ;
v_0 = 0.8 ;

y = zeros( (tot_ts+1), 1) ;
y(1) = 0.5 ;
y(2) = 0.5 + v_0*tau ;

for k = 1:tot_ts
    y(k+2) = (y(k)*(gamma*tau/2.0-1)+y(k+1)*(2.0-w*w*tau*tau))/(1+gamma*tau/2.0) ;
end

[t,yy]=ode45(@damped1,[0 30],[0.5 0.8]); 
figure,plot(t,yy(:,1)), hold on, plot(0:tau:(tot_time+tau),y,'r')
%figure,plot(t,yy(:,2)), hold on, plot(0:tau:tot_time,y(:,2),'r')