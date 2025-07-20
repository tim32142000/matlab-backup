% numerical method: RK4
% last modified: 2015.01.21

clear all
tic
tau = 0.0001 ;
tot_time = 1000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
b2_per_ts = int32( (2*pi/w_ext)/tau) ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 103.4 ;
Odiv = 2 ; % divide a circle to how many part
delta = 0 ; % initial phase of external B field
alpha = pi/2 ; % angle from B_1 direction to B_2 up-direction
%delta = rand(1,1)*2*pi ;
%delta_de = delta/pi*180 ;

theta = zeros( (tot_ts+1),Odiv) ;
w = zeros( (tot_ts+1),Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;
w_n = zeros(floor(tot_time),Odiv) ;
%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
wetau = w_ext*tau ;
for k = 1:2
    theta(1,k) = 0.0*2*pi + -1*(k-1)*(k-3)*0.2778*2*pi + (k-1)/2.0*(k-2)*0.66389*2*pi ;
    %theta(1,k) = 0.0 + -1*(k-1)*(k-3)*1.57 + (k-1)/2.0*(k-2)*0.66389*2*pi ;
    w(1,k) = 0.0 ;
    jj = 0 ; % indices of theta_n and <theta>_n
    for m = 1:tot_ts % use rk4 to solve ODEs
        rkdth1 = w(m,k) ;
        rkdw1 = -gamma*rkdth1 - b1*sin(theta(m,k) ) + b2*sin(alpha-theta(m,k) )*cos( m*wetau+delta) ; % have been corrected
        rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
        rkdw2 = -gamma*rkdth2 - b1*sin(theta(m,k) + 0.5*rkdth1*tau) + b2*sin(alpha-(theta(m,k) + 0.5*rkdth1*tau) )*cos( (m+0.5)*wetau+delta) ; % have been corrected
        rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
        rkdw3 = -gamma*rkdth3 - b1*sin(theta(m,k) + 0.5*rkdth2*tau) + b2*sin(alpha-(theta(m,k) + 0.5*rkdth2*tau) )*cos( (m+0.5)*wetau+delta) ; % have been corrected
        rkdth4 = rkdth1 + tau*rkdw3 ;
        rkdw4 = -gamma*rkdth4 - b1*sin(theta(m,k) + rkdth3*tau) + b2*sin(alpha-(theta(m,k) + rkdth3*tau) )*cos( (m+1)*wetau+delta) ; % have been corrected
        theta(m+1,k) = theta(m,k) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
        w(m+1,k) = w(m,k) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
        if(mod(m,b2_per_ts)==0)
            jj = jj + 1 ;
            theta_n(jj,k) = theta(m+1,k) ;
            w_n(jj,k) = w(m+1,k) ;
            nc = floor((theta_n(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_n(jj,k) = theta_n(jj,k) - nc*2*pi ;
            end
            theta_ave(jj,k) = mean(theta( (m+1-b2_per_ts+1):m+1,k) ) ;
            nc = floor((theta_ave(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj,k) = theta_ave(jj,k) - nc*2*pi ;
            end
        end
    end
    for m = 1:tot_ts
        nc = floor( (theta(m+1,k) + pi)/2/pi) ;
        if( nc ~= 0)
                theta(m+1,k) = theta(m+1,k) - nc*2*pi ;
        end
    end

   
    figure; plot(theta( ( (tot_time-300)/tau):10:(tot_time/tau),k)./(2*pi),w( (tot_time-300)/tau:10:(tot_time/tau),k)./(2*pi),'.','MarkerSize',2)
    xlabel('\theta(x2\pi)')
    ylabel('\omega(x2\pi/T)')
    %ylabel('$\dot \theta(\times2\pi/T)$','interpreter','latex','fontsize',20)
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
   
    
    %{
    hold on
    plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    hold off
    
    figure; plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', num2str(b2),', \theta_0=',num2str(theta(k,1)/pi*180 ),'\circ'])
    %}
    
    %{
    figure; plot(theta_n(tot_time-300:tot_time-1,k)./(2*pi),theta_n(tot_time-299:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('\theta_n')
    ylabel('\theta_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure; plot(theta_ave(tot_time-300:tot_time-1,k)./(2*pi),theta_ave(tot_time-299:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('<\theta>_n')
    ylabel('<\theta>_{n+1}')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure;plot( 0:tau:20,theta(tau/tau:20/tau+1,k)./(2*pi),'b')
    xlabel('time(T)')
    ylabel('\theta')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])

    figure;plot( (tot_time-20):tau:tot_time,theta( (tot_time-20)/tau:tot_time/tau,k)./(2*pi),'b')
    xlabel('time(T)')
    ylabel('\theta')
    title(['B_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/pi*180.0 )])
    ylim([-0.5 0.5])
    %}
    
end

%{
b0=sqrt(b1*b1 + b2*b2);
for k=1:16
    time=(tot_time-10+(k-1)/16.0);
    figure
    compass(exp(1i*theta(time/tau,1)),'r')
    hold on
    compass(b1/b0,b2*cos(2*pi*time)/b0,'k');
    compass(exp(1i*theta(time/tau,2)),'b')
    hold off
    title(['B_2=' num2str(b2,'%.2f') ' t=' num2str(time,'%.1f') '(T)'])
end
%}
runtime = toc