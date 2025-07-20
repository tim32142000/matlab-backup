% Calculating Basin of Attraction(BA) by RK4 method
% last modified: 2015.03.30

clear all
tic

tau = 0.0001 ;
tot_time_max = 1000 ;
tot_ts_max = round(tot_time_max/tau) ; % ts:time step
w_ext = 2*pi ;
b2_per_ts = round(2*pi/w_ext/tau) ; % per:period
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 103.4 ;
Odiv = 100 ; % divide a circle to how many part
w_min_de = -620 ;
w_max_de = 620 ;
w_intv_de = 10 ;

%baplot = zeros(3,(floor( (w_max-w_min)/w_intv) + 1)*Odiv ) ; % 1st row which kind attractor, 2nd row theta_0, 3rd row w_0
nbaplota = 0 ; % count of 1st kind period attractor
nbaplotb = 0 ; % count of 2nd kind period attractor
coss = 10 ; % condition of successive satisfying
cssa = 0 ; % current successive satisfying for a
cssb = 0 ; % current successive satisfying for b

theta_n1a = -0.356629306087072 ; % a condition. from other program, in radius
theta_n1b = -1.259731006067219 ; % b condition

wetau = w_ext*tau ;
nthw = 0 ; % index of initial conditions
    
%figure; hold on;
for p = w_min_de:w_intv_de:w_max_de % loop change w_0
    toc
    fprintf('Start w=%3.0f\n',p)
    for k = 1:Odiv % loop change theta_0
        %nthw = nthw + 1 ;
        fprintf('Start Odiv=%3.0f of %3.0f\n',k,Odiv)
        
        clear w theta theta_n
        w = zeros(1, (tot_ts_max+1) ) ;
        theta = zeros(1, (tot_ts_max+1) ) ;
        theta_n = zeros(1,tot_time_max) ;
        theta_ave = zeros(1,tot_time_max) ;
        
        w(1) = p/180*pi ;
        theta(1) = -pi + (k-1)/Odiv*2*pi ;          
        %baplot(2,nthw) = theta(1) ;
        %baplot(3,nthw) = w(1) ;
        
        jj = 0 ; % count of theta_n
        m = 0 ; % count of while loop
        leavewhilem = 0 ; % switch variable: 0:had not satisfied stable condition; 1:satisfied 
        while (m <= tot_ts_max) && (leavewhilem == 0)
            m = m+1 ;
            rkdth1 = w(m) ;
            rkdw1 = -gamma*w(m) - b1*sin(theta(m) ) + b2*cos(theta(m) )*cos( m*wetau) ; % have been corrected
            rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
            rkdw2 = -gamma*rkdth2 - b1*sin(theta(m) + 0.5*rkdth1*tau) + b2*cos(theta(m) + 0.5*rkdth1*tau)*cos( (m+0.5)*wetau) ; % have been corrected
            rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
            rkdw3 = -gamma*rkdth3 - b1*sin(theta(m) + 0.5*rkdth2*tau) + b2*cos(theta(m) + 0.5*rkdth2*tau)*cos( (m+0.5)*wetau) ; % have been corrected
            rkdth4 = rkdth1 + tau*rkdw3 ;
            rkdw4 = -gamma*rkdth4 - b1*sin(theta(m) + rkdth3*tau) + b2*cos(theta(m) + rkdth3*tau)*cos( (m+1)*wetau) ; % have been corrected
            theta(m+1) = theta(m) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
            w(m+1) = w(m) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
            if(mod(m,b2_per_ts)==0 )
                jj = jj+1 ;
                theta_n(jj) = theta(m+1) ;
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                if( nc ~= 0)
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                theta_ave(jj) = mean(theta((m+1-b2_per_ts+1):m+1) ) ;
                nc = floor((theta_ave(jj) + pi)/2/pi) ;
                if( nc ~= 0)
                    theta_ave(jj) = theta_ave(jj) - nc*2*pi ;
                end
                if  (abs(theta_n(jj) - theta_n1a) < 10^-7)
                    cssa = cssa + 1 ;
                    if (cssa >= coss)
                        nbaplota = nbaplota + 1 ;
                        baplota(1,nbaplota) = theta(1) ;
                        baplota(2,nbaplota) = w(1) ;
                        leavewhilem = 1 ;
                    end
                    %plot(theta(1)/pi*180,w(1)/pi*180,'color','b','LineStyle','none','Marker','.','MarkerSize',6) % blue stands for 1st kind of periodic attractor
                elseif (abs(theta_n(jj) - theta_n1b) < 10^-7)
                    cssb = cssb + 1 ;
                    if (cssb >= coss)
                        nbaplotb = nbaplotb + 1 ;
                        baplotb(1,nbaplotb) = theta(1) ;
                        baplotb(2,nbaplotb) = w(1) ;
                        leavewhilem = 1 ;
                    end
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
figure; hold on
plot(baplota(1,:)/pi*180,baplota(2,:)/pi*180,'color','b','LineStyle','none','Marker','.','MarkerSize',6) ;
plot(baplotb(1,:)/pi*180,baplotb(2,:)/pi*180,'color','r','LineStyle','none','Marker','.','MarkerSize',6) ;
hold off
xlim([-180-360/Odiv 180]);ylim([(w_min_de-w_intv_de) (w_max_de+w_intv_de)])
xlabel('\theta_0(\circ)')
ylabel('\omega_0(\circ/T)')
title(['B_2=', sprintf('%.2f',b2), ' RK4'])

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