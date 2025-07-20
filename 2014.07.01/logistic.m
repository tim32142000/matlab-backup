clear all
tot_steps = 1000 ;
r_min = 2.4;
r_max = 4.0 ;
r_intv=0.001 ;
r_num = int16( (r_max-r_min)/r_intv + 1) ;
x = zeros(r_num,tot_steps+1) ;

x(:,1) = 0.49 ;
r = r_min ;
figure; hold on;
for k = 0:r_num-1
    for n = 1:tot_steps
        x(k+1,n+1) = r*x(k+1,n)*(1-x(k+1,n) ) ;
    end
    plot(ones(1,600).*r,x(k+1,401:1000),'.','MarkerSize',2)
    r = r + r_intv ;
end
hold off;
xlabel('r')
