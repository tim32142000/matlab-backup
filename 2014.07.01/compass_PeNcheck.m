p = 2 ; % check pwriod-p
k = 1 ;

figure; 
semilogy(p+1:tot_time,abs(theta_n(p+1:tot_time)-theta_n(1:tot_time-p) ) )
title(['Euler B_2=', num2str(B2), ' \theta_0=', num2str(theta(k,1)/pi*180.0), '\circ \omega_0=', num2str(w(k,1)/pi*180.0), '(\circ/\tau) \tau=10^{', num2str(log10(tau) ),'}'])
% better to add numerical method used to title
xlabel('n(index for \theta_n)')
ylabel(['log_{10}(|\theta_n-\theta_{n-',num2str(p),'}|)'])
