tic
clear all
tau = 0.0001 ;
tot_time = 1000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
B1 = 36.0 ;
B2 = 99.30 ;
Odiv = 1 ; % divide a circle to how many part

w = zeros(Odiv, (tot_ts+1) ) ;
theta = zeros(Odiv, (tot_ts+1) ) ;
theta_n = zeros(Odiv,floor(tot_time) ) ;
%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for k = 1:1
    theta(k,1) = 0.0/180.0*pi ;
    w(k,1) = 0.0 ;
    jj = 1 ;
    for m = 1:tot_ts % use rk4 for solving ODE
        rkdth1 = w(k,m) ;
        rkdw1 = -gamma*w(k,m) - B1*sin(theta(k,m) ) - B2*cos(theta(k,m) )*cos( m*wetau) ;
        rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
        rkdw2 = -gamma*rkdth2 - B1*sin(theta(k,m) + 0.5*rkdth1*tau) - B2*cos(theta(k,m) + 0.5*rkdth1*tau)*cos( (m+0.5)*wetau) ;
        rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
        rkdw3 = -gamma*rkdth3 - B1*sin(theta(k,m) + 0.5*rkdth2*tau) - B2*cos(theta(k,m) + 0.5*rkdth2*tau)*cos( (m+0.5)*wetau) ;
        rkdth4 = rkdth1 + tau*rkdw3 ;
        rkdw4 = -gamma*rkdth4 - B1*sin(theta(k,m) + rkdth3*tau) - B2*cos(theta(k,m) + rkdth3*tau)*cos( (m+1)*wetau) ;
        theta(k,m+1) = theta(k,m) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
        w(k,m+1) = w(k,m) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
        %theta(k,m+2) = (theta(k,m)*(gata2-1.0) + theta(k,m+1)*2.0 - tausq*(B1*sin(theta(k,m+1) )+B2*cos(theta(k,m+1) )*cos(m*wetau) ) )/(1+gata2) ;
        %w(k,m+1) = (theta(k,m+2) - theta(k,m) )/(2*tau) ;
        if(mod(m,tor_ext_ts)==0)
            theta_n(k,jj) = theta(k,m+1) ;
            if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                theta_n(jj) = theta_n(jj) - nc*2*pi ;
            end
            jj = jj+1 ;
        end
    end
    for m = 1:tot_ts
        if( floor((theta(k,m+1) + pi)/2/pi) ~= 0)
                nc = floor((theta(k,m+1) + pi)/2/pi) ;
                theta(k,m+1) = theta(k,m+1) - nc*2*pi ;
        end
    end

    figure; plot(theta(k,( (tot_time-500)/tau):(tot_time/tau) )./pi*180,w(k,(tot_time-500)/tau:(tot_time/tau) )./pi*180,'.','MarkerSize',2)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', num2str(B2),', \theta_0=',num2str(theta(k,1)/pi*180 ),'\circ'])
end
figure; plot(theta_n(tot_time-500:tot_time-1)./pi*180,theta_n(tot_time-499:tot_time)./pi*180,'.','MarkerSize',6)
xlabel('\theta_n')
ylabel('\theta_{n+1}')
title(['B_2=', num2str(B2),', \theta_0=',num2str(theta(k,1)/pi*180 ),'\circ'])

figure;plot( (tot_time-20):tau:tot_time,theta( (tot_time-20)/tau:tot_time/tau)./pi*180,'b')
xlabel('time')
ylabel('\theta(\circ)')
%{
figure; plot(theta(1,(200/tau):(401/tau)-1),w(1,200/tau:(401/tau)-1),...
    theta(2,(200/tau):(401/tau)-1),w(2,200/tau:(401/tau)-1),...
    theta(3,(200/tau):(401/tau)-1),w(3,200/tau:(401/tau)-1),...
    theta(4,(200/tau):(401/tau)-1),w(4,200/tau:(401/tau)-1),...
    theta(5,(200/tau):(401/tau)-1),w(5,200/tau:(401/tau)-1),...
    theta(6,(200/tau):(401/tau)-1),w(6,200/tau:(401/tau)-1),...
    theta(7,(200/tau):(401/tau)-1),w(7,200/tau:(401/tau)-1),...
    theta(8,(200/tau):(401/tau)-1),w(8,200/tau:(401/tau)-1),...
    'LineStyle','none','Marker','.','MarkerSize',4)
axis([-pi pi -15 15])
%}
%{
%[t,yy]=ode45(@compass_odefn,[0 101],[0.5 1],[],B1,B2); 
figure%,plot(t,yy(:,1)), hold on
plot(0:tau:(tot_time+tau),theta,'b')
%figure,plot(t,yy(:,2)), hold on, plot(0:tau:tot_time,y(:,2),'r')
figure; plot(theta((80/tau):(601/tau)-1),w(80/tau:(601/tau)-1),'.')
%plot(yy(7000:,1),yy(7000:,2) );
%}
runtime = toc