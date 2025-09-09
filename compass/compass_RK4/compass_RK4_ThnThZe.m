% Plot Theta_N versus Theta_0

clear all
tic
tau = 0.0001 ;
tot_time = 300 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
delta = 0 ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 95.0 ;
Odiv = 2 ; % divide a circle to how many part

%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for n = 1:1
    b2=b2 ;
    %figure; hold on;
    for k = 1:Odiv
        clear w theta theta_n
        w = zeros(1, (tot_ts+1)) ;
        theta = zeros(1, (tot_ts+1)) ;
        theta_n = zeros(1,tot_time) ;
        w(1) = 0.0 ;
        theta(1) = -pi + (k-1)/Odiv*2*pi ;
        jj = 1 ;
        for m = 1:tot_ts
            rkdth1 = w(m) ;
            rkdw1 = -gamma*w(m) - b1*sin(theta(m) ) - b2*cos(theta(m) )*cos( m*wetau+delta) ;
            rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
            rkdw2 = -gamma*rkdth2 - b1*sin(theta(m) + 0.5*rkdth1*tau) - b2*cos(theta(m) + 0.5*rkdth1*tau)*cos( (m+0.5)*wetau+delta) ;
            rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
            rkdw3 = -gamma*rkdth3 - b1*sin(theta(m) + 0.5*rkdth2*tau) - b2*cos(theta(m) + 0.5*rkdth2*tau)*cos( (m+0.5)*wetau+delta) ;
            rkdth4 = rkdth1 + tau*rkdw3 ;
            rkdw4 = -gamma*rkdth4 - b1*sin(theta(m) + rkdth3*tau) - b2*cos(theta(m) + rkdth3*tau)*cos( (m+1)*wetau+delta) ;
            theta(m+1) = theta(m) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
            w(m+1) = w(m) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
            if(mod(m,tor_ext_ts)==0 )
                theta_n(jj) = theta(m+1) ;
                if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                    nc = floor((theta_n(jj) + pi)/2/pi) ;
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                jj = jj+1 ;
            end
        end
    plot(theta(1)/pi*180*ones(1,300),theta_n(tot_time-299:tot_time)./pi.*180,'LineStyle','none','Marker','.','MarkerSize',6 )
    %figure; plot(theta(k,(500/tau):(700/tau)-1),w(k,500/tau:(700/tau)-1),'.')
    end
    %hold off;
    xlim([-180 180]);ylim([-180 180])
    xlabel('\theta_0(\circ)')
    ylabel('\theta_n(\circ)')
    title(['B_2=', num2str(b2)])
end
toc
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