phi = [ -16.3
-12
-7.6
-3.7
-1.3
0
0.8
2.8
6.7
10.4
15.6];
fre = [5700
5900
6100
6300
6400
6460
6510
6600
6800
7000
7300 ];

y = phi;
x = fre ;
% dx = (0.0001/sqrt(12))*ones(15,1);
% dy = (1/sqrt(12))*ones(15,1);
% dy_trasf = dx.*(3809);

inc_tot = d_fase_1(18:28);

    
w = 1./(inc_tot).^2;
Delta = ((sum(w))*(sum(w.*(x.^2))))-((sum(w.*x)).^2);
A = (((sum(w.*(x.^2)))*(sum(w.*y)))-((sum(w.*x))*(sum(w.*(x.*y)))))/(Delta);
B = (((sum(w))*(sum(w.*(x.*y))))-((sum(w.*y))*(sum(w.*x))))/(Delta);
sigma_A_quadro = (sum(w.*(x.^2)))/(Delta);
sigma_B_quadro = (sum(w))/(Delta);
sigma_A = sqrt(sigma_A_quadro);
sigma_B = sqrt(sigma_B_quadro);

chi2 = sum(((y-(A+B*x)).^2)./(inc_tot.^2))

hfig = fit_ab_graph(x, y, zeros(size(x)), inc_tot, A, B, 'Hz', 'phi(deg)', 'phi(deg)', 0, 1, 0, 1);


% model_x = B*x + A
% x_max = x(15)
% 


