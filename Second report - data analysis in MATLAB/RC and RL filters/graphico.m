R=[20.21 60.614 110.6 161.04 200.54];
R_inc=[0.206 0.21 0.21 0.22 0.224];
fine_dati = 6500;  
dati_cas = 300;
V_in = ones(dati_cas,1);
V_out = ones(dati_cas,1);
time = ones(dati_cas,1);

%% IMPORT DATI
for i= 1:25
filename = strcat('./dati/scope_',int2str((i-1)),'.csv');
Q{i}=csvread(filename,3,0,[3,0,9000,2]);
end

%% FACCIO UN CICLO PER OGNI RESISTENZA 

V_out_m = double(vpa(Q{1}(50:8500,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{1}(40:8500,1)));
min_time = min(time_m);

V_in = double(vpa(Q{1}(50:7000,2)))-min_V_out;
V_out= double(vpa(Q{1}(50:7000,3)))-min_V_out;
time = double(vpa(Q{1}(50:7000,1)))-min_time;

d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(1);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
alfa=fit_out(1);
beta=fit_out(2);
ceta=fit_out(3);

scatter(time,V_in,.9,'m');
grid on;
hold on;
scatter(time,V_out,1.1,'y');
scatter(time, exp(alfa+beta.*time+ceta./V_out),.9,'k');
plot(3.9e-6*[1 1],[-2.60 1.17],'m','LineWidth', 1);
legend('V_{in}','V_{out}', 'Regressione');
title('V_{in} e V_{out}','FontSize',13, 'FontName', 'David Libre');
yl = ylabel('V');
set(yl, 'FontSize', 9);
xl = xlabel('s');
set(xl, 'FontSize', 9);


M = [0.1:0.1:100];
lyap_m = [1:1:1000];
fig1=figure();
scatter(M,lyap_m, 1.5);   
title('Esponente di Lyapunov in funzione di \mu (\delta(0) =1e-8)','FontSize',13, 'FontName', 'David Libre');
yl = ylabel('\lambda');
set(yl, 'FontSize', 13);
xl = xlabel('\mu = m/M');
set(xl, 'FontSize', 13);
grid on
hold off






