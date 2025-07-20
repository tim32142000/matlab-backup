clear all
tic
tau = 0.0001 ;
tot_time = 1000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
B1 = 36.0 ;
B2 = 99.30 ;
Odiv = 1 ; % divide a circle to how many part

w = zeros(Odiv, (tot_ts+2) ) ;
theta = zeros(Odiv, (tot_ts+2) ) ;
theta_n = zeros(Odiv,floor(tot_time) ) ;
%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for k = 1:1
    w(k,1) = 0.0 ;
    w(k,2) = 0.0 ;
    theta(k,1) = 288.0/180.0*pi ;
    theta(k,2) = theta(k,1) + w(k,1)*tau ;
    jj = 1 ;
    for m = 1:tot_ts
        theta(k,m+2) = (theta(k,m)*(gata2-1.0) + theta(k,m+1)*2.0 - tausq*(B1*sin(theta(k,m+1) )+B2*cos(theta(k,m+1) )*cos(m*wetau) ) )/(1+gata2) ;
        %{
        if( floor((theta(k,m+2) + pi)/2.0/pi) ~= 0)
            nc = floor( (theta(k,m+2) + pi)/2/pi) ;
            theta(k,m+2) = theta(k,m+2) - nc*2*pi ;
        end
        %}
        w(k,m+1) = (theta(k,m+2) - theta(k,m) )/(2*tau) ;
        if(mod(m,tor_ext_ts)==0)
            theta_n(k,jj) = theta(k,m+2) ;
            if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                theta_n(jj) = theta_n(jj) - nc*2*pi ;
            end
            jj = jj+1 ;
        end
    end
    for m = 1:tot_ts
        if( floor((theta(k,m+2) + pi)/2/pi) ~= 0)
                nc = floor((theta(k,m+2) + pi)/2/pi) ;
                theta(k,m+2) = theta(k,m+2) - nc*2*pi ;
        end
    end

    w(k,tot_ts+2) = (theta(tot_ts+2)-theta(tot_ts+1) )/tau ;
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