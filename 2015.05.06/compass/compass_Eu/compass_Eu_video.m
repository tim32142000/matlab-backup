% numerical method: Euler
% last modified: 2014.12.30
% unit of theta: radius/(2*pi)

clear all
tic
tau = 0.0001 ;
tot_time = 500 ;
tot_ts = uint32(tot_time/tau) ;
w_ext = 2*pi ;
b2_per_ts = uint32( (2*pi/w_ext)/tau) ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 100.6 ;
Odiv = 1 ; % divide a circle to how many part
nf = uint16(1) ; % index of figure

delta = 0.0 ; % initial phase of external B field
%delta = rand(1,1)*2*pi ;
%delta_de = delta/pi*180 ;

theta = zeros((tot_ts+2),Odiv) ;
w = zeros( 2,Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;

%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
if ( Odiv >= 1 )
    w(1,1) = 0.0 ;
    w(2,1) = 0.0 ;
    theta(1,1) = 0.25*2*pi ;
    theta(2,1) = theta(1,1) + w(1,1)*tau ;
end
if (Odiv >= 2 )
    w(1,2) = 0.0 ;
    w(2,2) = 0.0 ;
    theta(1,2) = 0.278*2*pi ;
    theta(2,2) = theta(1,2) + w(1,2)*tau ;
end
if (Odiv >= 3 )
    w(1,3) = 0.0 ;
    w(2,3) = 0.0 ;
    theta(1,3) = -0.3361*2*pi ;
    theta(2,3) = theta(1,3) + w(1,3)*tau ;
end

figure(nf)
switch Odiv
    case {1}
        for m = 1:double(tot_ts)
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_per_ts/50)==0 && m>2000000)
                compass(1.5*exp(1i*theta(m+2,1)),'r')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                hold off
                title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                %title(['real time:',num2str(toc,'%.2f'),'(sec) simulation time:',num2str(m*tau,'%.2f'),'(T)']);
                
                figure(nf)
            end

        end
    case {2}
        for m = 1:double(tot_ts)
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_per_ts/50)==0)
                compass(1.5*exp(1i*theta(m+2,1)),'b')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                compass(1.5*exp(1i*theta(m+2,2)),'r')
                hold off
                title(['real time:',num2str(toc,'%.2f'),'(sec) simulation time:',num2str(m*tau,'%.2f'),'(T)']);
                legend(['\theta_0=' num2str(theta(1,1)/(2*pi),'%.3f')],'B field',['\theta_0=' num2str(theta(1,2)/(2*pi),'%.3f')],'Location','NorthEastOutside')
                figure(nf)
                %pause(0.1)
            end

        end
    case {3}
        for m = 1:double(tot_ts)
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_per_ts/50)==0)
                compass(1.5*exp(1i*theta(m+2,1)),'b')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                compass(1.5*exp(1i*theta(m+2,2)),'r')
                compass(1.5*exp(1i*theta(m+2,3)),'g')
                hold off
                title(['real time:',num2str(toc,'%.2f'),'(sec) simulation time:',num2str(m*tau,'%.2f'),'(T)']);
                figure(nf)
            end
        end
end
hold off
nf = nf + 1 ;




runtime = toc