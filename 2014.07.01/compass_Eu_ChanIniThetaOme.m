clear all
tic
tau = 0.0001 ;
tot_time_max = 400 ;
tot_ts_max = tot_time_max/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
B1 = 36.0 ;
B2 = 99.3 ;
Odiv = 10 ; % divide a circle to how many part
w_min = -600.0/180.0*pi ;
w_max = 600.0/180.0*pi ;
w_intv = 50.0/180.0*pi ;

theta_n2a = [1.395779873278010 1.228883903717313] ; % from other program, in radius
theta_n2b = [0.505826174094475 0.050935426265787] ;

%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for n = 1:1 % loop change B2
    B2=B2 ;
    figure; hold on;
    for p = -600:50:600
        toc
        fprintf('Start w=%3.0f\n',p)
        for k = 1:Odiv %loop change theta_0
            fprintf('Start Odiv=%3.0f of %3.0f\n',k,Odiv)
            clear w theta theta_n
            w = zeros(1, 2) ;
            theta = zeros(1, (tot_ts_max+1)) ;
            theta_n = zeros(1,tot_time_max) ;
            w(1) = p/180*pi ;
            w(2) = 0.0 ;
            theta(1) = (k-1)/Odiv*2*pi ;
            theta(2) = theta(1) + w(1)*tau ;
            jj = 1 ;
            m = 1 ; % count of while loop
            leavewhilem = 0 ;
            while m <= tot_ts_max && leavewhilem == 0
                theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 - tausq*(B1*sin(theta(m+1) )+B2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
                if(mod(m,tor_ext_ts)==0 )
                    theta_n(jj) = theta(m+2) ;
                    if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                        nc = floor((theta_n(jj) + pi)/2/pi) ;
                        theta_n(jj) = theta_n(jj) - nc*2*pi ;
                    end
                    if ( (abs(theta_n(jj) - theta_n2a(2) ) < 10^-6) && (abs(theta_n(jj-1) - theta_n2a(1) ) < 10^-6) )
                        leavewhilem = 1 ;
                        plot(theta(1)/pi*180,w(1)/pi*180,'color','b','LineStyle','none','Marker','.','MarkerSize',6) % blue stands for 1st kind of periodic attractor
                    elseif ( (abs(theta_n(jj) - theta_n2b(2) ) < 10^-6) & (abs(theta_n(jj-1) - theta_n2b(1) ) < 10^-6) )
                        leavewhilem = 2 ;
                        plot(theta(1)/pi*180,w(1)/pi*180,'color','r','LineStyle','none','Marker','.','MarkerSize',6) % red stands for 2nd kind of periodic attractor
                    end
                    jj = jj + 1 ;
                end
                m = m + 1 ;
            end
            
    %figure; plot(theta(k,(500/tau):(700/tau)-1),w(k,500/tau:(700/tau)-1),'.')
        end
    end
    hold off;
    xlim([0 360]);ylim([-700 700])
    xlabel('\theta_0(\circ)')
    ylabel('\omega_0(\circ/\tau)')
    title(['B_2=', num2str(B2)])
end
runtime = toc
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