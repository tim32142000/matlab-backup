% numerical method: Euler
% last modified: 2015.04.24

clear all
tic
tau = 0.0001 ;
tot_time = 1000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
b2_per_ts = int32( (2*pi/w_ext)/tau) ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 102.8 ;
Odiv = 2 ; % divide a circle to how many part
delta = 0 ; % initial phase of external B field
%delta = rand(1,1)*2*pi ;
%delta_de = delta/pi*180 ;

theta = zeros((tot_ts+2),Odiv) ;
w = zeros( (tot_ts+2),Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;
w_n = zeros(floor(tot_time),Odiv) ;

%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for k = 1:2
    w(1,k) = 4.945652205825546 + (k-1)*eps  ;
    w(2,k) = 0.0 ;
    %theta(1,k) = 0.0*2*pi + -1*(k-1)*(k-3)*0.25*2*pi ;
    theta(1,k) = -1.222045690793823 + -1*(k-1)*(k-3)*eps ;
    theta(2,k) = theta(1,k) + w(1,k)*tau ;
    jj = 1 ; % indices of theta_n and <theta>_n
    for m = 1:tot_ts
        theta(m+2,k) = (theta(m,k)*(gata2-1.0) + theta(m+1,k)*2.0 + tausq*(-b1*sin(theta(m+1,k) )+b2*sin(pi/2-0.01-theta(m+1,k) )*cos(m*wetau+delta) ) )/(1+gata2) ;
        %{
        if( floor((theta(k,m+2) + pi)/2.0/pi) ~= 0)
            nc = floor( (theta(k,m+2) + pi)/2/pi) ;
            theta(k,m+2) = theta(k,m+2) - nc*2*pi ;
        end
        %}
        w(m+1,k) = (theta(m+2,k) - theta(m,k) )/(2*tau) ;
        if(mod(m,b2_per_ts)==0)
            theta_n(jj,k) = theta(m+2,k) ;
            w_n(jj,k) = w(m+1,k) ;
            nc = floor((theta_n(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_n(jj,k) = theta_n(jj,k) - nc*2*pi ;
            end
            theta_ave(jj,k) = mean(theta( (m+2-b2_per_ts+1):m+2,k) ) ;
            nc = floor((theta_ave(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj,k) = theta_ave(jj,k) - nc*2*pi ;
            end
            jj = jj+1 ;
        end
    end
    for m = 1:tot_ts
        nc = floor((theta(m+2,k) + pi)/2/pi) ;
        if( nc ~= 0)
            theta(m+2,k) = theta(m+2,k) - nc*2*pi ;
        end
    end

    w(tot_ts+2,k) = (theta(tot_ts+2,k)-theta(tot_ts+1,k) )/tau ;
    
   
    
    figure; plot(theta( ( (tot_time-400)/tau):5:(tot_time/tau),k)./(2*pi),w( (tot_time-400)/tau:5:(tot_time/tau),k)./(2*pi),'.','MarkerSize',2)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    
    
    %{
    hold on
    plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    hold off
    
    figure; plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', num2str(b2),', \theta_0=',num2str(theta(k,1)/pi*180 ),'\circ'])
    %}



    figure; plot(theta_n(tot_time-200:tot_time-1,k)./(2*pi),theta_n(tot_time-199:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('\theta_n')
    ylabel('\theta_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure; plot(theta_ave(tot_time-200:tot_time-1,k)./(2*pi),theta_ave(tot_time-199:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('<\theta>_n')
    ylabel('<\theta>_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure;plot( 0:tau:20,theta(tau/tau:20/tau+1,k)./(2*pi),'b')
    xlabel('time')
    ylabel('\theta')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])

    figure;plot( (tot_time-20):tau:tot_time,theta( (tot_time-20)/tau:tot_time/tau,k)./(2*pi),'b')
    xlabel('time')
    ylabel('\theta')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])

end

runtime = toc