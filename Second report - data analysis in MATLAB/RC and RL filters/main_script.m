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
Q{i}=csvread(filename,300,0,[300,0,9700,2]);
end

%% FACCIO UN CICLO PER OGNI RESISTENZA 
for i=1:5
numeri = int64(rand([1 dati_cas])*fine_dati);
V_out_m = double(vpa(Q{i}(1:end,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{i}(1:end,1)));
min_time = min(time_m);

for h=1:dati_cas
V_in(h) = double(vpa(Q{i}(numeri(h)+1,2)));
V_out(h)= double(vpa(Q{i}(numeri(h)+1,3)))-min_V_out;
time(h) = double(vpa(Q{i}(numeri(h)+1,1)))-min_time;
end
d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(3);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
H_1{i}=-fit_out(2);
end

%%
for i=6:10
numeri = int64(rand([1 dati_cas])*fine_dati);
V_out_m = double(vpa(Q{i}(1:end,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{i}(1:end,1)));
min_time = min(time_m);

for h=1:dati_cas
V_in(h) = double(vpa(Q{i}(numeri(h)+1,2)));
V_out(h)= double(vpa(Q{i}(numeri(h)+1,3)))-min_V_out;
time(h) = double(vpa(Q{i}(numeri(h)+1,1)))-min_time;
end
d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(3);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
H_2{i-5}=-fit_out(2);
end
%%
for i=11:15
numeri = int64(rand([1 dati_cas])*fine_dati);
V_out_m = double(vpa(Q{i}(1:end,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{i}(1:end,1)));
min_time = min(time_m);

for h=1:dati_cas
V_in(h) = double(vpa(Q{i}(numeri(h)+1,2)));
V_out(h)= double(vpa(Q{i}(numeri(h)+1,3)))-min_V_out;
time(h) = double(vpa(Q{i}(numeri(h)+1,1)))-min_time;
end
d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(3);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
H_3{i-10}=-fit_out(2);
end
%%
for i=16:20
numeri = int64(rand([1 dati_cas])*fine_dati);
V_out_m = double(vpa(Q{i}(1:end,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{i}(1:end,1)));
min_time = min(time_m);

for h=1:dati_cas
V_in(h) = double(vpa(Q{i}(numeri(h)+1,2)));
V_out(h)= double(vpa(Q{i}(numeri(h)+1,3)))-min_V_out;
time(h) = double(vpa(Q{i}(numeri(h)+1,1)))-min_time;
end
d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(3);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
H_4{i-15}=-fit_out(2);
end
%%
for i=21:25
numeri = int64(rand([1 dati_cas])*fine_dati);
V_out_m = double(vpa(Q{i}(1:end,3)));
min_V_out=min(V_out_m);
time_m = double(vpa(Q{i}(1:end,1)));
min_time = min(time_m);

for h=1:dati_cas
V_in(h) = double(vpa(Q{i}(numeri(h)+1,2)));
V_out(h)= double(vpa(Q{i}(numeri(h)+1,3)))-min_V_out;
time(h) = double(vpa(Q{i}(numeri(h)+1,1)))-min_time;
end
d_logV = log(ones(size(V_out))*0.5*8*3/100/2);
d_time = ones(size(V_out))*(8e-04)*4.5*0.005/R(3);
[fit_out, dfit_out, C, chi2, N_DOF] = lsq_fit_gen(log(V_out),[ones(size(V_out)) time 1./V_out],'err', d_logV);
H_5{i-20}=-fit_out(2);
end

%% IN QUESTA SEZIONE IL PEDICE INDICA IL VALORE DI RESISTENZA 
B_1=mean(cell2mat(H_1));            
B_2=mean(cell2mat(H_2));
B_3=mean(cell2mat(H_3));
B_4=mean(cell2mat(H_4));
B_5=mean(cell2mat(H_5));
B_tot=[B_1 B_2 B_3 B_4 B_5];

std_1=std(cell2mat(H_1))/sqrt(5);
std_2=std(cell2mat(H_2))/sqrt(5);
std_3=std(cell2mat(H_3))/sqrt(5);
std_4=std(cell2mat(H_4))/sqrt(5);
std_5=std(cell2mat(H_5))/sqrt(5);
inc_tot=[std_1 std_2 std_3 std_4 std_5];

%% FIT FINALE CON TEST CHI QUADRATO
y=B_tot;
x=R;
inc_tot_1 = sqrt(inc_tot.^2+(120*R_inc).^2)
w = 1./(inc_tot_1).^2;
Delta = ((sum(w))*(sum(w.*(x.^2))))-((sum(w.*x)).^2);
A = (((sum(w.*(x.^2)))*(sum(w.*y)))-((sum(w.*x))*(sum(w.*(x.*y)))))/(Delta)
B = (((sum(w))*(sum(w.*(x.*y))))-((sum(w.*y))*(sum(w.*x))))/(Delta)
sigma_A_quadro = (sum(w.*(x.^2)))/(Delta);
sigma_B_quadro = (sum(w))/(Delta);
sigma_A = sqrt(sigma_A_quadro)
sigma_B = sqrt(sigma_B_quadro)

chi2_1 = sum(((y-(A+B*x)).^2)./(inc_tot_1.^2))

    
    
    
    
    