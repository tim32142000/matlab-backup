clear all
tot_steps = 1000 ;
a_min = 0.9;
a_max = 4.0 ;
a_intv=0.001 ;
a_num = int16( (a_max-a_min)/a_intv + 1) ;
x = zeros(a_num,tot_steps+1) ;

x(:,1) = 0.49 ;
a = a_min ;
figure; hold on;
for k = 0:a_num-1
    for n = 1:tot_steps
        x(k+1,n+1) = x(k+1,n)*(a*x(k+1,n)*x(k+1,n)+(1-a) ) ;
    end
    plot(ones(1,600).*a,x(k+1,401:1000),'.','MarkerSize',2)
    a = a + a_intv ;
end
hold off;
xlabel('r')
