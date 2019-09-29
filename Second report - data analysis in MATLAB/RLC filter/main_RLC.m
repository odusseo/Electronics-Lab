R=[500.64 10002];
R_osc = 1002350;
R_l = 12;
C = 58.01e-9;
C_osc = 105.6e-12;
L = 0.01036;
C_l = (1/(4*pi*pi*132000*132000*L))-C_osc;
fine_dati = 40;

%% IMPORT DATI
for i= 1:2
filename = strcat('./res_',int2str(i),'.csv');
Q{i}=csvread(filename,2,0,[2,0,41,5]);
end

%% definire proprietà filtro RC, frequenze di studio (valori usati in aula)

f_1 = double(vpa(Q{1}(1:fine_dati,1))); 
f_2 = double(vpa(Q{2}(1:fine_dati,1))); 

Vin_1 = double(vpa(Q{1}(1:fine_dati,2))); 
Vin_2 = double(vpa(Q{2}(1:fine_dati,2))); 

Vout_1 = double(vpa(Q{1}(1:fine_dati,3))); 
Vout_2 = double(vpa(Q{2}(1:fine_dati,3))); 

fase_deg_1 = double(vpa(Q{1}(1:fine_dati,4))); 
fase_deg_2 = double(vpa(Q{2}(1:fine_dati,4))); 

fase_rad_1 = -deg2rad(fase_deg_1);
fase_rad_2 = -deg2rad(fase_deg_2);

dt_1 = 8e-4*10/2/1000*double(vpa(Q{1}(1:fine_dati,5))); %8e-4*10 sez orizz / 2 sigma / 1000 per trasf ms-->s* dtempo per sezione 
dt_2 = 8e-4*10/2/1000*double(vpa(Q{2}(1:fine_dati,5)));

d_fase_1 = 2*dt_1*2*pi.*f_1;
d_fase_2 = 2*dt_2*2*pi.*f_2;

dv_in_12= 3*8/2/100/1000*500;
dv_out_1 = 3*8/2/100/1000*double(vpa(Q{1}(1:fine_dati,6))); %(3.../100) full scale *8 sez vert/ 2 sigma / 1000 per trasf mV-->V* dtempo per sezione 
dv_out_2 = 3*8/2/100/1000*double(vpa(Q{2}(1:fine_dati,6)));

dv_H_1 = sqrt((dv_out_1./Vin_1).^2+(dv_in_12*(Vout_1)./(Vin_1.^2)).^2);
dv_H_2 = sqrt((dv_out_2./Vin_2).^2+(dv_in_12*(Vout_2)./(Vin_2.^2)).^2);

omega_1= 2*pi*f_1; 
omega_2= 2*pi*f_2; 

H_exp_1 = Vout_1./Vin_1;
H_exp_2 = Vout_2./Vin_2;

C_tot= C+C_osc+C_l;
f_3dB = 1./(2*pi*sqrt(L*(C+C_osc+C_l)))
dC_tot = sqrt((0.4e-9)^2+(1.6e-12)^2+(6.8e-12)^2);
df_3db_teo = sqrt((0.00003/(4*pi*sqrt((L^3)*C_tot)))^2+(dC_tot/(4*pi*sqrt(L*(C_tot)^3)))^2)

%% calcolare funzione di trasferimento
f = 10.^[1:.0001:6]';
omega = 2*pi*f; 

Z_c = 1./(j*omega*(C+C_osc+C_l));
Z_l = R_l+j*omega*L;
Z_par = 1./(1./Z_c + 1/R_osc + 1./Z_l);

H_1 = Z_par./ (R(1)+Z_par);
H_2 = Z_par./ (R(2)+Z_par);

Z_cp = 1./(j*omega*(C_osc+C_l));
Z_lp = R_l+j*omega*L;
Z_parp = 1./( 1/R_osc + 1./Z_lp);

H_1p = Z_parp./ (R(1)+Z_parp);
H_2p = Z_parp./ (R(2)+Z_parp);

H_prova_1 =(R_l/(R(1)+R_l))*ones(50001,1);
H_prova_2 =1./(1+j*omega*R(1)*(C+C_osc+C_l));


%% RES 1
fig_num1 = bode_plot(f,H_1,'add_dB');
axes(fig_num1(2));
plot(f_3dB*[1 1],[0.001 1],'k','LineWidth', 2)
axes(fig_num1(3));
plot(f_3dB*[1 1],[90 -90],'k', 'LineWidth', 2)
fig_num1 = bode_plot(f_1,[H_exp_1,fase_rad_1],'add_dB', 'points','.', ... 
    'err',[dv_H_1,d_fase_1],'fig',fig_num1,'ylim',[5e-3 2], 'col', 'r', 'Ohm');
hfig1 = gcf;
legend('Modello teorico','f_0', 'Dati sperimentali');
title(fig_num1(2), 'Bode plot - R_p','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');


%% RES 2
fig_num2 = bode_plot(f,H_2,'add_dB');
axes(fig_num2(2));
plot(f_3dB*[1 1],[0.001 1],'k','LineWidth', 2)
axes(fig_num2(3));
plot(f_3dB*[1 1],[90 -90],'k','LineWidth', 2)
fig_num2 = bode_plot(f_2,[H_exp_2, fase_rad_2],'add_dB', 'points','.', ... 
    'err',[dv_H_2, d_fase_2],'fig',fig_num2,'ylim',[1e-3 2], 'col', 'r');
hfig2 = gcf;
legend('Modello teorico','f_0', 'Dati sperimentali');
title(fig_num2(2), 'Bode plot - R_g','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');

%% RES 1 PROVA RL
fig_num1 = bode_plot(f,H_1p,'add_dB');
axes(fig_num1(2));
plot(f_3dB*[1 1],[0.001 1],'k','LineWidth', 2)
axes(fig_num1(3));
plot(f_3dB*[1 1],[90 -90],'k', 'LineWidth', 2)
fig_num1 = bode_plot(f_1,[H_exp_1,fase_rad_1],'add_dB', 'points','.', ... 
    'err',[dv_H_1,d_fase_1],'fig',fig_num1,'ylim',[5e-3 2], 'col', 'r');
hfig1 = gcf;
legend('Modello teorico','f_0', 'Dati sperimentali');
title(fig_num1(2), 'Bode plot - R_p - LR passa alto','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');

%% RES 1 PROVA RC
fig_num2 = bode_plot(f,H_prova_2,'add_dB');
axes(fig_num2(2));
plot(f_3dB*[1 1],[0.001 1],'k','LineWidth', 2)
axes(fig_num2(3));
plot(f_3dB*[1 1],[90 -90],'k', 'LineWidth', 2)
fig_num2 = bode_plot(f_1,[H_exp_1,fase_rad_1],'add_dB', 'points','.', ... 
    'err',[dv_H_1,d_fase_1],'fig',fig_num2,'ylim',[5e-3 2], 'col', 'r');
hfig2 = gcf;
legend('Modello teorico','f_0', 'Dati sperimentali');
title(fig_num2(2), 'Bode plot - R_p - RC passa basso','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');
yl = ylabel('V');
set(yl, 'FontSize', 9);
xl = xlabel('\Omega');
set(xl, 'FontSize', 9);




%%EXPORT figure

% Set paper orientation
set(hfig1, 'PaperOrientation', 'landscape');

% Set paper type
set(hfig1, 'PaperType', 'A4');

% Set paper positioning units
set(hfig1, 'PaperUnits', 'centimeters');

% Set paper positioning to be manual
set(hfig1, 'PaperPositionMode', 'manual');

% Set paper positioning to be filling better the page
set(hfig1, 'PaperPosition', [0 0 29 21]);

% exports a figure to a PNG format file, myfigure06.png
print(hfig1, '-dpng', 'myfigure01.png')

%%%%
% Set paper orientation
set(hfig2, 'PaperOrientation', 'landscape');

% Set paper type
set(hfig2, 'PaperType', 'A4');

% Set paper positioning units
set(hfig2, 'PaperUnits', 'centimeters');

% Set paper positioning to be manual
set(hfig2, 'PaperPositionMode', 'manual');

% Set paper positioning to be filling better the page
set(hfig2, 'PaperPosition', [0 0 29 21]);

% exports a figure to a PNG format file,
print(hfig2, '-dpng', 'myfigure02.png')


