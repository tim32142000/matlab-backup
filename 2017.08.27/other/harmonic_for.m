clear all
tau = 0.001 ;
tot_time = 100 ;
tot_ts = tot_time/tau ;
w = 2*pi ;
y = zeros(2,(tot_ts+1) ) ;

y(1,1) = 1 ;
y(2,1) = 1.5 ;

for k = 1:tot_ts
    y(2,k+1) = y(2,k) - tau*w*w*y(1,k) ;
    y(1,k+1) = y(1,k) + tau*y(2,k+1) ;
end

figure; plot(0:tau:tot_time, y(1,:) ) ; hold on ; plot(0:tau:tot_time, y(2,:) , 'r')
