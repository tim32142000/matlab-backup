% Calculating Basin of Attraction(BA) by Euler method

clear all
tic
tau = 0.0001 ;
tot_time_max = 500 ;
tot_ts_max = tot_time_max/tau ;
w_ext = 2*pi ;
b2_per_ts = int32(2*pi/w_ext/tau) ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 99.3 ;
Odiv = 360 ; % divide a circle to how many part
w_min_de = -650 ;
w_max_de = 650 ;
w_intv_de = 10 ;

w_min = w_min_de/180.0*pi ;
w_max = w_max_de/180.0*pi ;
w_intv = w_intv_de/180.0*pi ;

%baplot = zeros(3,(floor( (w_max-w_min)/w_intv) + 1)*Odiv ) ; % 1st row which kind attractor, 2nd row theta_0, 3rd row w_0
nbaplota = 0 ; % number of 1st kind period attractor
nbaplotb = 0 ; % number of 2nd kind period attractor
coss = 10 ; % condition of successive satisfying
cssa = 0 ; % current successive satisfying for a
cssb = 0 ; % current successive satisfying for b

theta_n2a = [1.395779873278010 1.228883903717313] ; % from other program, in radius
theta_n2b = [0.505826174094475 0.050935426265787] ;

%{
theta(1) = 0.69*pi ; % -0.5 and 0.5
theta(2) = theta(1) + w(1)*tau ;
%}
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;

%nthw = 0 ; % index of initial conditions
    
%figure; hold on;
for p = w_min_de:w_intv_de:w_max_de % loop change w_0
    toc
    fprintf('Start w=%3.0f\n',p)
    for k = 1:Odiv % loop change theta_0
        %nthw = nthw + 1 ;
        fprintf('Start Odiv=%3.0f of %3.0f\n',k,Odiv)
            
        clear w theta theta_n theta_ave
        w = zeros(1, 2) ;
        theta = zeros(1, (tot_ts_max+1)) ;
        theta_n = zeros(1,tot_time_max) ;
        theta_ave = zeros(1,tot_time_max) ;
            
        w(1) = p/180.0*pi ;
        %w(2)
        theta(1) = -pi + (k-1)/Odiv*2*pi ;
        theta(2) = theta(1) + w(1)*tau ;            
        %baplot(2,nthw) = theta(1) ;
        %baplot(3,nthw) = w(1) ;
           
        jj = 0 ; % count of theta_n
        m = 0 ; % count of while loop
        leavewhilem = 0 ; % switch variable: 0:had not satisfied stable condition; 1:had satisfied
        while m <= tot_ts_max && leavewhilem == 0
            m = m+1 ;
            theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 + tausq*(-b1*sin(theta(m+1) )+b2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
            if(mod(m,b2_per_ts)==0 )
                jj = jj+1 ;
                theta_n(jj) = theta(m+2) ;
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                if( nc ~= 0)
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                if ( (abs(theta_n(jj) - theta_n2a(2) ) < 10^-6) && (abs(theta_n(jj-1) - theta_n2a(1) ) < 10^-6) )
                    cssa = cssa + 1 ;
                    if (cssa >= coss)
                        nbaplota = nbaplota + 1 ;
                        baplota(1,nbaplota) = theta(1) ;
                        baplota(2,nbaplota) = w(1) ;
                        leavewhilem = 1 ;
                    end
                elseif ( (abs(theta_n(jj) - theta_n2b(2) ) < 10^-6) && (abs(theta_n(jj-1) - theta_n2b(1) ) < 10^-6) )
                    cssb = cssb + 1 ;
                    if (cssb >= coss)
                        nbaplotb = nbaplotb + 1 ;
                        baplotb(1,nbaplotb) = theta(1) ;
                        baplotb(2,nbaplotb) = w(1) ;
                        leavewhilem = 1 ;
                    end
                    %baplot(1,nthw) = 2 ;
                    %plot(theta(1)/pi*180,w(1)/pi*180,'color','r','LineStyle','none','Marker','.','MarkerSize',6) % red stands for 2nd kind of periodic attractor
                else
                    cssa = 0 ;
                    cssb = 0 ;
                end
            end
        end
            
    %figure; plot(theta(k,(500/tau):(700/tau)-1),w(k,500/tau:(700/tau)-1),'.')
    end
end
    
%hold off;
%{
nbaplota = 0 ;
nbaplotb = 0 ;

for m = 1:nthw % plot loop
    if baplot(1,m) == 1
        nbaplota = nbaplota + 1 ;
        baplota(1,nbaplota) = baplot(2,m) ;
        baplota(2,nbaplota) = baplot(3,m) ;
    elseif baplot(1,m) == 2
        nbaplotb = nbaplotb + 1 ;
        baplotb(1,nbaplotb) = baplot(2,m) ;
        baplotb(2,nbaplotb) = baplot(3,m) ;
    end
end
%}
figure; hold on
plot(baplota(1,:)/pi*180,baplota(2,:)/pi*180,'color','b','LineStyle','none','Marker','.','MarkerSize',6) ;
plot(baplotb(1,:)/pi*180,baplotb(2,:)/pi*180,'color','r','LineStyle','none','Marker','.','MarkerSize',6) ;
hold off
xlim([-180-360/Odiv 180]);ylim([(w_min_de-w_intv_de) (w_max_de+w_intv_de)])
xlabel('\theta_0(\circ)')
ylabel('\omega_0(\circ/time)')
title(['B_2=', sprintf('%.2f',b2),' Euler'])
        


runtime = toc
