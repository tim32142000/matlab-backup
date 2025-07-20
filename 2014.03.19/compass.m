clear all
tau = 0.0001 ;
tot_time = 701 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
ts_T_ext = 2*pi/w_ext/tau ;
gamma = 6 ;
B1 = 36 ;
B2 = 104 ;
Odiv = 8 ; % divide a circle to how many part

w = zeros(Odiv, (tot_ts+1)) ;
theta = zeros(Odiv, (tot_ts+1)) ;
%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for k = 1:Odiv
    w(k,1) = 0.0 ;
    w(k,2) = 0.0 ;
    theta(k,1) = (k-1)*2.0*pi/Odiv-pi ;
    theta(k,2) = theta(k,1) + w(k,1)*tau ;
    for l = 1:tot_ts
        theta(k,l+2) = (theta(k,l)*(gata2-1.0) + theta(k,l+1)*2.0 - tausq*(B1*sin(theta(k,l+1) )+B2*cos(theta(k,l+1) )*cos(l*wetau) ) )/(1+gata2) ;
        w(k,l+1) = (theta(k,l+2) - theta(k,l) )/(2*tau) ;
    end
    w(k,tot_ts+2) = (theta(tot_ts+2)-theta(tot_ts+1) )/tau ;
    figure; plot(theta(k,(3000/tau):(4000/tau)-1),w(k,3000/tau:(4000/tau)-1),'.')
end

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