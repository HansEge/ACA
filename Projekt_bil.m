close all;
clear;

load('bil_data.mat'); %Brug denne!

%Få data ud af gemte filer
motor_Deg = Rotation_Deg.signals.values;
motor_velocity = Velocity_Rads.signals.values;
stepInput = input.signals.values;

%Moteren vender forkert, så input er nødt til at være -100, derfor tages
%abs af værdierne.
ts = 0.01;
stepInput = abs(stepInput);
motor_velocity = abs(motor_velocity);
motor_Deg = abs(motor_Deg);

sys_deg = iddata(motor_Deg,stepInput,'Ts',ts);
sys_velocity = iddata(motor_velocity,stepInput,'Ts',ts);

tf_degrees = tfest(sys_deg,2,0);
tf_velocity = tfest(sys_velocity,1,0);

figure
plot(motor_Deg);
title('Motor deg');

figure
plot(motor_velocity);
title('Motor velocity');

[A,B,C,D]=tf2ss(tf_degrees.Numerator,tf_degrees.Denominator);

sysSS=ss(A,B,C,D);
%% Define requirements and desired location of poles
OS = 5;
Ts = 2;
zeta = -log(OS/100)/sqrt(pi^2+(log(OS/100)^2)); %desired damping ratio
wn = 4/zeta/Ts; %desired natural frequency

s = tf('s');
ch_eqn = wn^2/(s^2+2*zeta*wn*s+wn^2); %Charactaristic equation
poles = pole(ch_eqn);
K = place(sysSS.A,sysSS.B,poles);
%% Define closed loop system
A_cl = [A-B*K];
sysSS_cl = ss(A_cl,B,C,D);

figure
step(sysSS_cl);
title('Step response of sysSS_c_l');
%% Insert 3rd pole to eliminate ss error
p3 = -30;
poles_new = [poles' p3];

A_new = [A [0;0];-C 0];
B_new = [B;0];

K_new = place(A_new,B_new,poles_new);
Ke = -K_new(3);
K_new = [K_new(1) K_new(2)];

A_cl_new = [A-B*K_new B*Ke;-C 0];
B_cl_new = [0;0;1];
C_cl_new = [C 0];

sysSS_new = ss(A_cl_new,B_cl_new,C_cl_new,D);

figure
step(sysSS_new);
title('Step response of sysSS_n_e_w with 3rd pole');
%% Oberserver design in discrete form
sysSSD=c2d(sysSS,ts,'zoh');%Discrete

ch_eqnD = c2d(ch_eqn,ts,'tustin');
polesD = pole(ch_eqnD);

K_z = place(sysSSD.A,sysSSD.B,polesD);
A_cl = [sysSSD.A-sysSSD.B*K_z];

sysssD_cl = ss(A_cl,sysSSD.B,sysSSD.C,sysSSD.D,ts);

observer_req = poles*30;
observer_reqZ = polesD/30;

L = (place(sysSS.A',sysSS.C',observer_req))';
Ld = (place(sysSSD.A',sysSSD.C',observer_reqZ))';
%% Eliminate ss error in z-domain

p3_z = -0.000425;%HVORFOR ER SKAL DENNE POL LIGGE PRÆCIS I -0.000425???
polesD_new = [polesD' p3_z];

A_newD = [sysSSD.A [0;0];-sysSSD.C 0];
B_newD = [sysSSD.B;0];

K_newD = place(A_newD,B_newD,polesD_new);
Ke_z = -K_newD(3);

K_newD = [K_newD(1) K_newD(2)];

A_cl_newD = [sysSSD.A-sysSSD.B*K_newD sysSSD.B*Ke_z;-sysSSD.C 0];
B_cl_newD = [0;0;1];
C_cl_newD = [sysSSD.C 0];

sysss_newD = ss(A_cl_newD,B_cl_newD,C_cl_newD,D,ts);
figure
step(sysss_newD);
title('sysSS_n_e_w_D with Ke in z-domain');