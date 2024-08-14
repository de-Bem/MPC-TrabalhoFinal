Ts = 0.04;

N1 = 1;
N2 = 40; % horizonte de predição
Nuv = N2; % horizonte de controle
Nuw = N2;
e = 0.001; % incremento da derivação numerica para cálculo de G

pos_inicial = [0 0 0];

load('track.mat');
Referencia_x = Track_x;
Referencia_y = Track_y;
Referencia_theta = Track_theta;

delta_x = 1; % peso de referencia relativo a x
delta_y = 1;
delta_theta = 0;

lambda_v = 0.1; % peso de controle relativo a v
lambda_w = 0.01;

psi_x = 10; % peso da variável de folga (epsilon) relativo a x
psi_y = 10;
psi_theta = 10;

simTim = 10; % [s]
tol = 0.2; % tolerancia para fim da simulação

% RESTRIÇÕES:
% --------------------
x_max = 2;
y_max = 2.5;
theta_max = 5*pi;

x_min = -1;
y_min = 0;
theta_min = -theta_max;
% --------------------
v_max = 0.5; % m/s
w_max = pi/2; % rad/s

v_min = -v_max;
w_min = -w_max;
% --------------------
dv_max = Ts*v_max*10; % pode alcançar velocidade maxima em 1/10 segundo
dw_max = Ts*w_max*10;

dv_min = -dv_max;
dw_min = -dw_max;