Ts = 0.05;

N1 = 1;
N2 = 40; % horizonte de predição
Nuv = N2; % horizonte de controle
Nuw = N2;

e = 0.001; % incremento da derivação numerica para cálculo de G
stdev = 0.388; % desvio padrão do erro de modelagem (0.388 gera uma normal com 99% dos valores no intervalo [-1,1])

load('track.mat');

pos_inicial = [Track_x(1),Track_y(1),Track_theta(1)];

Referencia_x = Track_x;
Referencia_y = Track_y;
Referencia_theta = Track_theta;

delta_x = 1; % peso de referencia relativo a x
delta_y = 1;
delta_theta = 0;

lambda_v = 0.01; % peso de controle relativo a v
lambda_w = 0.01;

psi_x = 10; % peso da variável de folga (epsilon) relativo a x
psi_y = 10;
psi_theta = 10;

simTim = 20; % [s]
tol = 0.2; % tolerancia para fim da simulação

% RESTRIÇÕES:
% --------------------
x_max = max(Track_x);
y_max = max(Track_y);
theta_max = 5*pi;

x_min = min(Track_x);
y_min = min(Track_y);
theta_min = -theta_max;
% --------------------
v_max = 0.5; % m/s
w_max = pi/2; % rad/s

v_min = -v_max;
w_min = -w_max;
% --------------------
dv_max = Ts*v_max*10; % pode alcançar velocidade maxima em 1/15 segundos
dw_max = Ts*w_max*10;

dv_min = -dv_max;
dw_min = -dw_max;