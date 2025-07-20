clear all

tau = 0.001 ;
tot_time = 800 ;
tot_ts = tot_time/tau ;
w = 1 ;
b = 1.05 ;
g_min = 0 ;
g_max = 1 ;
n_g = 1 ;
g=linspace(g_min,g_max,n_g) ;

n_ini_phi = 20 ;
n_ini_theta = 20 ;
phi = zeros( (tot_ts+2),n_ini_phi) ;
theta = zeros( (tot_ts+2),n_ini_theta) ;


phi(1,:) = 0 ;
theta(1,:) = 0 ;

n_fig = gcf ;
figure(n_fig+1)
hold on

for k=1:n_g
    for m=1:tot_ts
        phi(m+1,k) = phi(m,k) + (1-b*mod(phi(m,k),1)+g*sin((theta(m,k)-phi(m,k))) )*tau;
        theta(m+1,k) = theta(m,k) + 1*tau ;
    end
    
    figure(n_fig+1)
    plot( (tot_time-100:100*tau:tot_time),phi((tot_time-100)/tau:100:(tot_time/tau),k) )
end
hold off

figure
