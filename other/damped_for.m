clear all
tau = 0.0001 ;
tot_time = 30 ;
tot_ts = tot_time/tau ;
w = 2*pi ;
gamma = 0.01 ;
y = zeros( (tot_ts+1), 2) ;

y(1,1) = 0.5 ;
y(1,2) = 0.8 ;

for k = 1:tot_ts
    y(k+1, 2) = y(k, 2) - tau*(gamma*y(k, 2) + w*w*y(k, 1) ) ;
    y(k+1, 1) = y(k, 1) + tau*y(k+1, 2) ;
end

[t,yy]=ode45(@damped1,[0 30],[0.5 0.8]); 
figure,plot(t,yy(:,1)), hold on, plot(0:tau:tot_time,y(:,1),'r')
%figure,plot(t,yy(:,2)), hold on, plot(0:tau:tot_time,y(:,2),'r')