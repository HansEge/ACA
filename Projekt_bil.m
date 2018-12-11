close all;

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

sysss=ss(A,B,C,D);

%% Define requirements and desired location of poles
OS = 5;
Ts = 2;
zeta = -log(OS/100)/sqrt(pi^2+(log(OS/100)^2)); %desired damping ratio
wn = 4/zeta/Ts; %desired natural frequency

s = tf('s');
ch_eqn = s^2+2*zeta*wn*s+wn^2;
[p,z,gain] = zpkdata(ch_eqn);
poles = cell2mat(p)'; %desired 2nd order poles for the closed loop system
K = place(sysss.A,sysss.B,poles);




%% Define closed loop system
A_cl = [A-B*K];
sysss_cl = ss(A_cl,B,C,D);

figure
step(sysss_cl);
title('Step response of sysss_c_l without 3rd pole');

%% Insert 3rd pole to eliminate ss error
p3 = -30;
poles_new = [poles p3];

A_new = [A [0;0];-C 0];
B_new = [B;0];

K_new = place(A_new,B_new,poles_new);
Ke = -K_new(3);
K_new = [K_new(1) K_new(2)];

A_cl_new = [A-B*K_new B*Ke;-C 0];
B_cl_new = [0;0;1];
C_cl_new = [C 0];

sysss_new = ss(A_cl_new,B_cl_new,C_cl_new,D);
figure
step(sysss_new);
title('Step response of sysss_n_e_w with 3rd pole');

%% Oberserver design in descrete form

sysssD=c2d(sysss,ts,'tustin');%Discrete

ch_eqnz = c2d(ch_eqn,ts,'tustin');
[p_z,z_z,gain_z] = zpkdata(ch_eqnz);
poles_z = cell2mat(p_z)';

K_z = place(sysssD.A,sysssD.B,poles_z);
A_cl = [sysssD.A-sysssD.B*K_z];

sysssD_cl = ss(A_cl,sysssD.B,sysssD.C,sysssD.D,ts);

observer_req = poles*5;

L = (place(sysss.A',sysss.C',observer_req))';

