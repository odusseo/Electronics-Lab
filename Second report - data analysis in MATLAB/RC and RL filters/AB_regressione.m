y = phi;
x = fre ;
% dx = (0.0001/sqrt(12))*ones(15,1);
% dy = (1/sqrt(12))*ones(15,1);
% dy_trasf = dx.*(3809);

inc_tot = 2*8e-4*10/2*0.000007576*2*pi.*fre

    
w = 1./(inc_tot).^2;
Delta = ((sum(w))*(sum(w.*(x.^2))))-((sum(w.*x)).^2);
A = (((sum(w.*(x.^2)))*(sum(w.*y)))-((sum(w.*x))*(sum(w.*(x.*y)))))/(Delta);
B = (((sum(w))*(sum(w.*(x.*y))))-((sum(w.*y))*(sum(w.*x))))/(Delta);
sigma_A_quadro = (sum(w.*(x.^2)))/(Delta);
sigma_B_quadro = (sum(w))/(Delta);
sigma_A = sqrt(sigma_A_quadro);
sigma_B = sqrt(sigma_B_quadro);

chi2 = sum(((y-(A+B*x)).^2)./(inc_tot.^2))

% model_x = B*x + A
% x_max = x(15)
% 


