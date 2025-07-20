% numerical method: Euler

clear all
tic
tau = 0.0001 ;
tot_time = 2000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
b1 = 36.0 ;
b2_min = 90.00 ;
b2_max = 90.00 ; % 117.77 over one circle
b2_intv = 0.01 ;
b2_num = int32( (b2_max-b2_min)/b2_intv+1 ) ;


%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
b2 = b2_min - b2_intv ;
w(1) = 0.0 ;
w(2) = 0.0 ;
theta = zeros(tot_ts+2,1) ;
theta(1) = 85.0/180.0*pi ;
theta(2) = theta(1) + w(1)*tau ;
theta_ave = zeros(1,tot_time) ;
thltavepl = zeros(2,b2_num) ; %theta long time average plot

%{
figure; 
hold on;
%}
for n = 1:b2_num
    b2 = b2 + b2_intv ;
    fprintf('Start b2=%6.2f, %5.0f of %5.0f\n',b2,n,b2_num)
    toc
    thltavepl(1,n) = b2 ;
    for m = 1:tot_ts
        theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 - tausq*(b1*sin(theta(m+1) )-b2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
    end
    thltavepl(2,n) = mean(theta( (tot_time/2)/tau:tot_time/tau) ) ;
    if( floor( (thltavepl(2,n) + pi)/(2*pi) ) ~= 0)
        nc = floor( (thltavepl(2,n) + pi) /(2*pi) ) ;
        thltavepl(2,n) = thltavepl(2,n) - nc*2*pi ;
    end
    %plot(b2*ones(1,500),theta_ave(501:1000)./pi.*180,'LineStyle','none','Marker','.','MarkerEdgeColor','r','MarkerSize',2 )
    %figure; plot(theta(k,(500/tau):(700/tau)-1),w(k,500/tau:(700/tau)-1),'.')
    %{
    w(1) = 0.0 ;
    w(2) = 0.0 ;
    theta(1) = 0.0 ;
    theta(2) = theta(1) + w(1)*tau ;
    %}
    
    w(1) = w(length(w)-1) ;
    w(2) = w(length(w)) ;
    theta(1) = theta(length(theta)-1) ;
    theta(2) = theta(length(theta)) ;
    
end
figure; plot(thltavepl(1,:),thltavepl(2,:)./pi*180.0,'LineStyle','none','Marker','.','MarkerEdgeColor','r','MarkerSize',6)
xlim([b2_min-b2_intv b2_max+b2_intv])
toc
%{
hold off;
xlabel('B_2')
ylabel('\theta_n(\circ)')
%}
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