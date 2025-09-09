clear all
tot_steps = 2000 ;
k_min = 2.0 ;
k_max = 6.0 ;
k_intv = 2^-6 ;
k_num = floor( (k_max-k_min)/k_intv + 1) ;
theta = zeros(tot_steps+1,1) ;
thpl_num = floor(600) ; % theta plot number
thpl = zeros(thpl_num*k_num,2) ; % theta plot

%rng(33)
%theta(1,1) = rand(k_num,1) ;
%theta(1,1) = 0.5 ;
theta0 = [0.25 0.5 0.75 (sqrt(5)-1)/2 ] ;
for m = 1:1
    for n = 1:k_num
        k = k_min + (n-1)*k_intv ;
        theta(1) = theta0(m) ;
        for o = 1:tot_steps
            theta(o+1) = mod(theta(o)-k/(2*pi)*sin(2*pi*theta(o) ),1) ;
        end
        thpl( (n-1)*thpl_num+1:n*thpl_num,1) = k ;
        thpl( (n-1)*thpl_num+1:n*thpl_num,2) = theta(tot_steps+1-thpl_num+1:tot_steps+1) ;
    end
    figure; plot(thpl(:,1),thpl(:,2),'.','MarkerSize',2)
    xlabel('K')
end


