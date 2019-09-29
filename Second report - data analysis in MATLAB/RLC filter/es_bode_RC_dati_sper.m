% esempi uso bode_plot con dati sperimentali...
% 
% per chiamare "help" per informazione:
% help bode_plot


% esempio semplice RC
f_3dB = 1e3;
f = 10.^[1:.001:5]';
H_RC = 1 ./ (1 + j*f/f_3dB);
fig = bode_plot(f,H_RC,'add_dB','col','r');

%%
% esempio dati sperimentali:
% dati completamente finti, quest'è solo un esempio per illustrare l'uso di bode plot:

f_exp = 10.^[1.5:0.5:5]';
H_mod_exp = [          1.032
    1.145
    0.998
    0.791
    0.243
    0.132
    0.0384
    0.0100];

H_phase_exp_deg =   [ 4.8
    3.0
    -13.9
    -42.9
    -73.4
    -85.1
    -80.5
    -88]; % in gradi

% convertire in radianti
H_phase_exp = H_phase_exp_deg * pi / 180;
     
dH_mod_exp =[         0.199900074937555
         0.199007438041998
         0.190692517849119
          0.14142135623731
        0.0603022689155527
        0.0199007438041998
       0.00632139541241014
       0.00199990000749938];

   dH_phase_exp =[                       0.1
                       0.1
                       0.1
                       0.1
                       0.1
                       0.1
                       0.1
                       0.1] * pi/180;
                   
%%
fig = bode_plot(f_exp,[H_mod_exp,H_phase_exp],'add_dB', ... 
    'points','.','err',[dH_mod_exp,dH_phase_exp],'fig',fig,'ylim',[5e-3 2]);