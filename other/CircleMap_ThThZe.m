clear all
tot_steps = 2000 ;
k_min = 2.0 ;
k_max = 6.0 ;
k_intv = 2^-6 ;
k_num = int32( (k_max-k_min)/k_intv + 1) ;
Odiv = 360 ; % devide circle to how many part
theta = zeros(tot_steps+1,Odiv) ;
theta(1,:)=(0:1/Odiv:1-1/Odiv) ;
thpl_num = int32(600) ; % theta plot number


for n = 1:k_num
    k = k_min + double((n-1))*k_intv ;
    
end