Ts = 0.05;

N1 = 1;
N2 = 30; % horizonte de predição
Nuv = 10; % horizonte de controle
Nuw = 25;
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

simTim = 45; % [s]
tol = 0.1; % tolerancia para fim da simulação

% RESTRIÇÕES:
% --------------------
x_max = 3;
y_max = 3;
theta_max = 10*pi;

x_min = -2;
y_min = -1;
theta_min = -theta_max;
% --------------------
v_max = 0.5;
w_max = pi/2;

v_min = -v_max;
w_min = -w_max;
% --------------------
dv_max = v_max;
dw_max = w_max;

dv_min = -dv_max;
dw_min = -dw_max;