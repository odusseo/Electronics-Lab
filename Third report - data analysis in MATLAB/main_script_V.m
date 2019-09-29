%% IMPORT DATI
clear all
close all
R_Z = [500
    290
    170
    80
    50
    30
    18.8
    16.8
    13.8
    10.5
    10.5
    10.5
    8.6
    8.6
    6.6
    6.6
    6.6
    5.6
    5.6
    5.6
    4.3
    4.3
    3.7
    3.7
    3.7
    ];

%I-V
filename = './diodo_rad_1.csv';
Q{1}=csvread(filename,2,0,[2,0,23,2]);
for i= 2:4
    filename = strcat('./diodo_rad_',int2str(i),'.csv');
    Q{i}=csvread(filename,2,0,[2,0,7,2]);
end
filename = './diodo_zen_1.csv';
Q{5}=csvread(filename,2,0,[2,0,31,2]);

%PONTE
C = 220e-6;
dC = 20/100*C;%0.4*C/100+0.1*10e-6/100;
filename = './ponte/ponte_1.csv';
Q{6}=csvread(filename,1,0,[1,0,27,5]);
filename = './ponte/ponte_2.csv';
Q{7}=csvread(filename,1,0,[1,0,25,7]);

%V TRASFORMATORE
filename = './VINO/scope_0.csv';
VIN_0=csvread(filename,10,0,[10,0,7600,1]);
filename = './VINO/scope_1.csv';
VIN_1=csvread(filename,10,0,[10,0,7600,1]);
filename = './VINO/scope_2.csv';
VIN_2=csvread(filename,10,0,[10,0,7600,1]);
filename = './VINO/scope_3.csv';
VIN_3=csvread(filename,10,0,[10,0,7600,1]);
filename = './VINO/scope_4.csv';
VIN_4=csvread(filename,10,0,[10,0,7600,1]);
filename = './VINO/scope_5.csv';
VIN_5=csvread(filename,10,0,[10,0,7600,1]);

%% ASEGNAZIONE VARIABILI
%I-V
di_1 = double(vpa(Q{1}(1:end,1)))/50/sqrt(12)/1000;
di_2 = double(vpa(Q{2}(1:end,1)))/50/sqrt(12)/1000;
di_3 = double(vpa(Q{3}(1:end,1)))/50/sqrt(12)/1000;
di_4 = double(vpa(Q{4}(1:end,1)))/50/sqrt(12)/1000;
di_5 = double(vpa(Q{5}(1:end,1)))/50/sqrt(12)/1000;

i_1 = double(vpa(Q{1}(1:end,2)))/1000;
i_2 = double(vpa(Q{2}(1:end,2)))/1000;
i_3 = double(vpa(Q{3}(1:end,2)))/1000;
i_4 = double(vpa(Q{4}(1:end,2)))/1000;
i_5 = double(vpa(Q{5}(1:end,2)))/1000;

v_1 = double(vpa(Q{1}(1:end,3)))/1000;
v_2 = double(vpa(Q{2}(1:end,3)))/1000;
v_3 = double(vpa(Q{3}(1:end,3)))/1000;
v_4 = double(vpa(Q{4}(1:end,3)))/1000;
v_5 = double(vpa(Q{5}(1:end,3)))/1000;

dv_1 = v_1*0.005/100+0.0035/100;
dv_2 = v_2*0.005/100+0.0035/100;
dv_3 = v_3*0.005/100+0.0035/100;
dv_4 = v_4*0.005/100+0.0035/100;
dv_5 = v_5*0.0035/100+10*0.0007/100;

%PONTE SENZA ZENER
r_1 = double(vpa(Q{6}(1:end,1)));
rip_1 = double(vpa(Q{6}(1:end,2)));
vmax_1 = double(vpa(Q{6}(1:end,3)));
d_vmax_1 = double(vpa(Q{6}(1:end,4)))/100000*1.5*8;
d_rip_1 = double(vpa(Q{6}(1:end,5)))/100/1000*1.5*8;
media_1 = double(vpa(Q{6}(1:end,6)));
for i=23:27
    media_1(i) = vmax_1(i)-0.01;
end
d_media_1 =  double(vpa(Q{6}(1:end,4)))/100000*1.5*8;
for i=1:13
    d_r_1(i) = 0.01*r_1(i)/100+0.001*1000/100;
end
for i=14:21
    d_r_1(i) = 0.01*r_1(i)/100+0.001*10000/100;
end
for i=22:27
    d_r_1(i) = 0.01*r_1(i)/100+0.001*10000/100;
end

%PONTE CON ZENER
r_2 = double(vpa(Q{7}(1:end,1)));
rip_2 = double(vpa(Q{7}(1:end,2)));
rip_c = double(vpa(Q{7}(1:end,3)));
d_rip_2 = double(vpa(Q{7}(1:end,4)))/100000*1.5*8;
d_rip_c = double(vpa(Q{7}(1:end,5)))/100000*1.5*8;
vmax_2 = double(vpa(Q{7}(1:end,6)));
vmax_c = double(vpa(Q{7}(1:end,7)));
d_vmax_2 = double(vpa(Q{7}(1:end,8)))/100000*1.5*8;
d_vmax_c = double(vpa(Q{7}(1:end,8)))/100000*1.5*8;
for i=1:21
    d_r_2(i) = 0.01*r_2(i)/100+0.001*1000/100;
end
for i=22:23
    d_r_2(i) = 0.01*r_2(i)/100+0.001*10000/100;
end
for i=24:25
    d_r_2(i) = 0.01*r_2(i)/100+0.001*10000/100;
end


%VINO
VIN_MAX = [max(VIN_0(1:2000,2)) ...
    max(VIN_0(2000:4000,2)) max(VIN_0(4000:end,2)) max(VIN_1(1:2000,2)) ...
    max(VIN_1(2000:4000,2)) max(VIN_1(4000:end,2)) max(VIN_2(1:2000,2)) ...
    max(VIN_2(2000:4000,2)) max(VIN_2(4000:end,2)) max(VIN_3(1:2000,2)) ...
    max(VIN_3(2000:4000,2)) max(VIN_3(4000:end,2)) max(VIN_4(1:2000,2)) ...
    max(VIN_4(2000:4000,2)) max(VIN_4(4000:end,2)) max(VIN_5(1:2000,2)) ...
    max(VIN_5(2000:4000,2)) max(VIN_5(4000:end,2))];
d_vinmax= sqrt((2.8*1.5*8/100).^2 +(std(VIN_MAX)./sqrt(15)).^2);
vinmax = 10.704;%mean(VIN_MAX);


%% %% %% %% %% %% %% ANALISI %% %% %% %% %% %% %% %%
%% I_V Calcolo rz
for i=2:30
    R_z(i-1)=abs((v_5(i)-v_5(i-1))./(i_5(i)-i_5(i-1)));
end

%% I-V PLOTS
fig0=figure();
figh0 = errorbar(v_1, i_1,di_1, di_1,dv_1,dv_1, 'Marker','.', 'Color', 'b', 'markersize',10,'LineStyle', 'none');
grid on
set(figh0, 'Capsize', 0)
hold on
x=[0.1:0.01:10];
plot(x,-2*x./r_1(1)+vinmax/r_1(1), 'color', 'r', 'LineWidth', 2);
xlim([0.1 5.5])
ylim([0 0.31])
yl = ylabel('A');
set(yl, 'FontSize', 18);
xl = xlabel('V');
set(xl, 'FontSize', 18);
llegend = legend( 'Curva i-V','Load line');
llegend.FontSize = 22;
rect = [0.65, 0.75, .1, .1];
set(llegend, 'Position', rect)
title('Curva i-V diodo e linea di carico','FontSize',17,'FontName', 'David Libre','FontWeight', 'normal');

hold off


% fig1=figure();
% hAy1=axes;
% hAy1.YScale='log';
% hold all
% figh1 = errorbar(v_1, i_1,di_1, di_1,dv_1,dv_1, 'LineStyle', 'none', 'Marker','.', 'Color', 'b');
% grid on
% set(figh1, 'Capsize', 0)
% hold off
% 
% 
% fig2 = figure();
% hAy2=axes;
% hAy2.YScale='log';
% hold all
% figh2 = errorbar(v_2, i_2,di_2, di_2,dv_2,dv_2, 'LineStyle', 'none', 'Marker','.', 'Color', 'b');
% grid on
% set(figh2, 'Capsize', 0)
% hold off
% 
% fig3 = figure();
% hAy3=axes;
% hAy3.YScale='log';
% hold all
% figh3 = errorbar(v_3, i_3,di_3, di_3,dv_3,dv_3, 'LineStyle', 'none', 'Marker','.', 'Color', 'b');
% grid on
% set(figh3, 'Capsize', 0)
% hold off
% 
% fig4 = figure();
% hAy4=axes;
% hAy4.YScale='log';
% hold all
% figh4 = errorbar(v_4, i_4,di_4, di_4,dv_4,dv_4, 'LineStyle', 'none', 'Marker','.', 'Color', 'b');
% grid on
% set(figh4, 'Capsize', 0)
% hold off

% zener
%for i=1:25
R=100;
fig5=figure();
figh5 = errorbar(v_5, i_5,di_5, di_5,dv_5,dv_5,  'Marker','.', 'Color', 'b', 'markersize',10,'LineStyle', 'none');
grid on
set(figh5, 'Capsize', 0)
hold on
x=[0.1:0.01:10];
plot(x,vmax_c(10)./R.*(-x.*(r_2(10)+R)./(r_2(10).*vmax_c(10))+1), 'color', 'r', 'LineWidth', 2);
xlim([2 7.5])
ylim([0 0.160])
yl = ylabel('A');
set(yl, 'FontSize', 18);
xl = xlabel('V');
set(xl, 'FontSize', 18);
llegend = legend( 'Curva i-V','Load line');
llegend.FontSize = 22;
rect = [0.25, 0.75, .1, .1];
set(llegend, 'Position', rect)
title('Curva i-V diodo e linea di carico','FontSize',17,'FontName', 'David Libre','FontWeight', 'normal');
hold off
%end 
% fig6=figure();
% hAy6=axes;
% hAy6.YScale='log';
% hold all
% figh6 = errorbar(v_5, i_5,di_5, di_5,dv_5,dv_5, 'LineStyle', 'none',
% 'Marker','.', 'Color', 'b');
% grid on
% set(figh6, 'Capsize', 0)
% hold off

%% PONTE NO ZENER
f=50;

for i=1:27
    syms omt1;
    omega = 2*pi*f;
    omt2 = atan(-omega*r_1(i)*C);
    OMT2(i) =  omt2;
    eqn = sin(omt1+pi) == sin(omt2)*exp(-(omt1+pi-omt2)/(omega*r_1(i)*C));
    solx= vpasolve(eqn, omt1);
    omt1 = double(solx);
    OMT1(i) =  omt1;
    
    I_dc(i) = vmax_1(i)/1097; %vmax_1(14)/(2*pi*1097)*(1-cos(omt2))
    E_r(i) = double(vpa((pi+omt1-omt2)*I_dc(i)/(omega*C)));

    E_dc_sing(i) = (vmax_1(i)/(2*pi))*(sqrt(1+omega*omega*r_1(i)*r_1(i)*C*C)*(1-cos(omt2-omt1)));
    E_dc_dopp(i) = vmax_1(i)-E_r(i)/2;
    E_dc_rms(i) = vpa(E_r(i)/(2*sqrt(3)));
end
lollete = (OMT1 - OMT2)/omega;
periodo = (lollete').*2;
f= 1./periodo;
inc_ripteo_1= sqrt((d_vmax_1./(f.*r_1.*C)).^2+(d_r_1.*vmax_1./(f.*r_1.*C.*r_1)).^2+(dC.*vmax_1./(f.*r_1.*C.*C)).^2);
inc_ripteo_1= inc_ripteo_1(1:end,1);
rip_1_teo = vmax_1./(2.*f.*C.*r_1);



fig7=figure();
hAy7=axes;
hAy7.YScale='log';
set(gca,'xscale','log')
hold all
figh7 = errorbar(r_1,rip_1,d_rip_1,d_rip_1, d_r_1,d_r_1, 'LineStyle', 'none', 'Marker','.','markersize',20, 'Color', 'b');
hold all
figh8 = errorbar(r_1,rip_1_teo,inc_ripteo_1,inc_ripteo_1, d_r_1,d_r_1, 'LineStyle', 'none', 'Marker','.', 'markersize',20, 'Color', 'r');
grid on
llegend = legend( 'Dati sperimentali','Modello teorico');
llegend.FontSize = 22;
rect = [0.65, 0.75, .1, .1];
set(llegend, 'Position', rect)
title('V_{ripple} vs R_L - logarithmic scale','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');
set(figh7, 'Capsize', 0)
set(figh8, 'Capsize', 0)
yl = ylabel('V');
set(yl, 'FontSize', 18);
xl = xlabel('\Omega');
set(xl, 'FontSize', 18);
hold off
% zoompos = [0.30   0.35    0.19    0.35];
zoompos = [0.15   0.15    0.19    0.35];
zoomaxs = axes('position', zoompos);
h_zoom = errorbar(r_1,rip_1_teo,inc_ripteo_1,inc_ripteo_1, d_r_1,d_r_1);
set(h_zoom, 'linestyle', 'none')
set(h_zoom, 'marker', '.', 'markersize', 16);
set(h_zoom, 'markeredgecolor', 'r')
set(h_zoom, 'markerfacecolor', 'r')
set(h_zoom, 'color', 'r')
set(h_zoom, 'capsize', 0)
grid on
hold on
idx = 10;
xlim(r_1(idx) + 1.5*[-d_r_1(idx) d_r_1(idx)])
set(zoomaxs, 'ylim', rip_1_teo(idx) + 1.5*[-inc_ripteo_1(idx) inc_ripteo_1(idx)]);

% OMT1 = OMT1';
% OMT2 = OMT2';
% OMT1 = rad2deg(OMT1);
% OMT2 = rad2deg(OMT2)+180;
for i=1:27
    syms z;
    a =   2.066e-23;
    b =       60.67;
    c =    4.31e-09;
    d =       21.29;
    eqn = a*exp(b*z) + c*exp(d*z) == -z./r_1(i)+vinmax/r_1(i);
    solx= vpasolve(eqn, z);
    inter(i) = double(solx);
end
inter = inter';
%vinmax = 7.5*sqrt(2)
d_vinmax = d_vinmax*ones(27,1);


% % zoompos = [0.30   0.35    0.19    0.35];
% zoompos = [0.15   0.15    0.19    0.35];
% zoomaxs = axes('position', zoompos);
% h_zoom = errorbar(r_1,rip_1_teo,inc_ripteo_1,inc_ripteo_1, d_r_1,d_r_1);
% set(h_zoom, 'linestyle', 'none')
% set(h_zoom, 'marker', '.', 'markersize', 16);
% set(h_zoom, 'markeredgecolor', 'r')
% set(h_zoom, 'markerfacecolor', 'r')
% set(h_zoom, 'color', 'r')
% set(h_zoom, 'capsize', 0)
% grid on
% hold on
% idx = 10;
% xlim(r_1(idx) + 1.5*[-d_r_1(idx) d_r_1(idx)])
% set(zoomaxs, 'ylim', vinmax-2*inter(idx) + 1.5*[-d_vinmax(idx) d_vinmax(idx)]);



%% PONTE ZENER
R=100;
for i=1:25
   syms t   
   a =   1.164e-08  ;
   b =       3.684 ;
   c =  -1.326e-08 ;
   d =       3.657  ;
   eqn = a*exp(b*t) + c*exp(d*t)== vmax_c(i)./R.*(-t.*(r_2(i)+R)./(r_2(i).*vmax_c(i))+1);
   solx= vpasolve(eqn, t);
   interz(i) = double(solx);
end
interz=interz';
%R_Z=abs(1./(R_Z'));

%(R_Z+R+(R.*R_Z./r_2))
rip_2_teo = rip_c.*R_Z./(R_Z+R+(R.*R_Z./r_2));
d_rip_2_teo = sqrt((d_rip_c).^2+(rip_c.*d_r_2.*R_Z.*(R_Z.*R./r_2.^2)./(R_Z+R+(R.*R_Z./r_2)).^2).^2+...
    (d_rip_c.*R_Z./(R_Z+R+(R.*R_Z./r_2))).^2);
d_rip_2_teo = d_rip_2_teo(1:end,1);

vd_2_teo = r_2./(r_2+R).*vmax_c;
d_d_2_teo = sqrt((d_vmax_c.*r_2./(r_2+R)).^2+(d_r_2.*vmax_c.*( R./(r_2+R).^2)).^2);
d_d_2_teo = d_d_2_teo(1:end,1);
%% rout
vmax_2 = vmax_2 -rip_2./2;
i_L = vmax_2./r_2;
r_out = ((5.26-0.032/2)-vmax_2)./i_L;
a = (d_r_2.*(5.27-vmax_2)./vmax_2).^2
b = (d_vmax_2.*r_2.*5.27./(vmax_2.*vmax_2)).^2
c = (vmax_2(25)./i_L).^2
d_r_out = sqrt((d_r_2.*(5.27-vmax_2)./vmax_2).^2+...
    (d_vmax_2.*r_2.*5.27./(vmax_2.*vmax_2)).^2+...
    (d_vmax_2(25)./i_L).^2);
d_r_out = d_r_out(1:end,1);


fig8=figure();
% hAy8=axes;
% hAy8.YScale='log';
set(gca,'xscale','log')
hold all
figh13 = errorbar(r_2,r_out,d_r_out,d_r_out, d_r_2,d_r_2, 'LineStyle', 'none', 'Marker','.', 'Color', 'b','markersize',20);
hold all
grid on
%figh14 = errorbar(r_2,interz, d_d_2_teo, d_d_2_teo, d_r_2,d_r_2, 'LineStyle', 'none', 'Marker','.', 'Color', 'r', 'markersize',20);
%llegend = legend( 'Dati sperimentali','Modello teorico');
%llegend.FontSize = 22;
%rect = [0.65, 0.25, .1, .1];
%set(llegend, 'Position', rect)
%xlim([90 50000])
%ylim([4.2 5.4])
title('R_{out} vs R_L','FontSize',15,'FontName', 'David Libre','FontWeight', 'normal');
yl = ylabel('R_{out} [\Omega]');
set(yl, 'FontSize', 18);
xl = xlabel('R_{L} [\Omega]');
set(xl, 'FontSize', 18);
hold off
set(figh13, 'Capsize', 0)
%set(figh14, 'Capsize', 0)



vout = vmax_c.*(4.5./(4.5+100*(1+(4.5*r_2))));

%% EXPORT figure

hfig1 = fig0;

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
% % Set paper orientation
% set(hfig2, 'PaperOrientation', 'landscape');
% 
% % Set paper type
% set(hfig2, 'PaperType', 'A4');
% 
% % Set paper positioning units
% set(hfig2, 'PaperUnits', 'centimeters');
% 
% % Set paper positioning to be manual
% set(hfig2, 'PaperPositionMode', 'manual');
% 
% % Set paper positioning to be filling better the page
% set(hfig2, 'PaperPosition', [0 0 29 21]);
% 
% % exports a figure to a PNG format file,
% print(hfig2, '-dpng', 'myfigure02.png')
% 


