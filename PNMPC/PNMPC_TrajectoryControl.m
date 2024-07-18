% Modelo cinematico de um robo:
%     x'=v*cos(theta)
%     y'=v*sin(theta)
%     theta'=w
% 
% x, y e theta representam a posição global do robo, enquanto v e w são 
% sinais de controle para velocidade linear e angular do robo.
% 
% O modelo discretizado:
% 
% x(k+1) = x(k) + T*v(k)*cos(theta(k))
% y(k+1) = y(k) + T*v(k)*sin(theta(k))
% theta(k+1) = theta(k) + T*w(k)
% 
% -> u(k) = u(k-1) + du(k)
% 
%     x(k+1) = x(k) + T*(v(k-1)+dv(k))*cos(theta(k))
%     y(k+1) = y(k) + T*(v(k-1)+dv(k))*sin(theta(k))
%     theta(k+1) = theta(k) + T*(w(k-1)+dw(k))

clc
close all
clear all

%% ------------------------- CONFIGURAÇÕES ------------------------------


Ts = 0.1;

N1 = 1;
N2 = 20; % horizonte de predição
Nuv = N2; % horizonte de controle
Nuw = N2;
e = 0.001; % incremento da derivação numerica para cálculo de G

pos_inicial = [0 0 0];

Referencia_x = 0.5;
Referencia_y = 1;
Referencia_theta = 0;

Qe_x = 0.8; % peso de referencia relativo a x
Qe_y = 5;
Qe_theta = 0.8;

Qu_v = 0; % peso de controle relativo a v
Qu_w = 0;

psi_x = 5; % peso da variável de folga (epsilon) relativo a x
psi_y = 5;
psi_theta = 5;

simTim = 20; % [s]
tol = 0.05; % tolerancia para fim da simulação

% RESTRIÇÕES:
% --------------------
x_max = 1;
y_max = 1;
theta_max = 2*pi;

x_min = 0;
y_min = 0;
theta_min = -theta_max;
% --------------------
v_max = 0.5;
w_max = pi/4;

v_min = -v_max;
w_min = -w_max;
% --------------------
dv_max = v_max/5;
dw_max = w_max/5;

dv_min = -dv_max;
dw_min = -dw_max;


%% --------------------------- MONTAGEM --------------------------------


N = N2-N1+1;

Ref_x = Referencia_x*ones(N,1);
Ref_y = Referencia_y*ones(N,1);
Ref_theta = Referencia_theta*ones(N,1);
Ref = [Ref_x; Ref_y; Ref_theta]; % transformando em vetor para fazer operações como -2*G_total'*Qe*(Ref-f_total)


Qe_x = Qe_x*eye(N)/(x_max^2*N); 
Qe_y = Qe_y*eye(N)/(y_max^2*N);
Qe_theta = Qe_theta*eye(N)/(theta_max^2*N);
Qe = [   Qe_x   zeros(N)  zeros(N);
       zeros(N)   Qe_y    zeros(N);
       zeros(N) zeros(N) Qe_theta]; % transformando em matriz para fazer operações como 2*(G_total'*Qe*G_total+Qu)


Qu_v = Qu_v*eye(Nuv)/(v_max^2*Nuv);
Qu_w = Qu_w*eye(Nuw)/(w_max^2*Nuw);
Qu = [     Qu_v      zeros(Nuv,Nuw);
      zeros(Nuw,Nuv)    Qu_w]; % transformando em matriz para fazer operações como 2*(G_total'*Qe*G_total+Qu)


Gthetav = zeros(N2,Nuv); % o sinal de velocidade não influencia a variável de estado theta
Gthetaw=[]; % o sinal de velocidade angular tem relação linear com a a variável de estado theta
for i = 1:N2
    for j = 1:Nuw
        if j <= i
            Gthetaw(i, j) = Ts*(i - j + 1);
        end
    end
end


%% Restrições


% RESTRIÇÃO DE CONTROLE
A_restr_controle = [tril(ones(Nuv))  zeros(Nuv,Nuw); 
                   -tril(ones(Nuv))  zeros(Nuv,Nuw);
                    zeros(Nuw,Nuv)   tril(ones(Nuw)); 
                    zeros(Nuw,Nuv)  -tril(ones(Nuw))];

b_aux_restr_controle = [v_max*ones(Nuv,1);
                       -v_min*ones(Nuv,1);
                        w_max*ones(Nuw,1);
                       -w_min*ones(Nuw,1)];

% -------------------------------------------------------------------------

% RESTRIÇÃO DE INCREMENTO DE CONTROLE

A_restr_increm = [eye(Nuv) zeros(Nuv,Nuw);
                 -eye(Nuv) zeros(Nuv,Nuw); 
                  zeros(Nuw,Nuv)  eye(Nuw); 
                  zeros(Nuw,Nuv) -eye(Nuw)];

b_restr_increm = [dv_max*ones(Nuv,1); -dv_min*ones(Nuv,1);
                  dw_max*ones(Nuw,1); -dw_min*ones(Nuw,1)];

% -------------------------------------------------------------------------

% RESTRIÇÃO DE SAÍDA
b_aux_restr_saida = [x_max*ones(N,1);  y_max*ones(N,1);  theta_max*ones(N,1);
                    -x_min*ones(N,1); -y_min*ones(N,1); -theta_min*ones(N,1)];

% VARIÁVEL DE FOLGA

% os ultimos 6 termos da solução do quadprog serão:
% [ eps_x_min, eps_y_min, eps_theta_min, eps_x_max, eps_y_max, eps_theta_max ]
n_epsilons = 6;

% matriz Aw a partir de A
n_rows_Aw = length(b_aux_restr_controle)+length(b_restr_increm)+length(b_aux_restr_saida);
Aw = zeros(n_rows_Aw,Nuv+Nuw+n_epsilons);

TR = zeros(N * n_epsilons, n_epsilons); % top right (G;-G é o primeiro bloco da matriz A)
for i = 1:n_epsilons
    TR((i-1)*N+1:i*N, i) = -1; % matrix com os epsilons para somar com G
end

Aw(1:size(TR,1) , end-n_epsilons+1:end) = TR;

% matriz Hw a partir de Hqp
BR = 2*diag([psi_x psi_y psi_theta psi_x psi_y psi_theta]); % bottom right

Hw = zeros(Nuw+Nuv+n_epsilons);
Hw(end-n_epsilons+1:end,end-n_epsilons+1:end) = BR;

% vetor fw a partir de fqp
fw = [zeros(Nuw+Nuv+n_epsilons,1)];


%% --------------------------- SIMULAÇÃO ----------------------------------


ini_sim = N1+1;
fim_sim = ini_sim+round(simTim/Ts);

dUs = zeros(fim_sim,2);
Us  = zeros(fim_sim,2);
Ys  = zeros(fim_sim,3);
Is  = zeros(fim_sim,Nuv+Nuw+n_epsilons);
Epsilons = zeros(fim_sim,n_epsilons);

% -------------------------------------------------------------------------

% necessário para fazer simulação encerrar de acordo com a tolerancia 
% especificada, se eu usar Ys(k,:) diretamente, o compilador entende que 
% norm(Ref) > tol, já que Ys foi inicializado com zeros
position_k = [0 0 0];

% estabelecendo condições iniciais
Ys(1:ini_sim,:) = repmat(pos_inicial,ini_sim,1);

% inicializações
resposta_livre = zeros(N2,3);
Gxv=[0];Gyv=[0];

% -------------------------------------------------------------------------

k=ini_sim;
tic
while (k<=fim_sim && norm([Referencia_x Referencia_y Referencia_theta]-position_k(1:3)) > tol)
    % calculando o modelo
    Ys(k,1) = Ys(k-1,1) + Ts*Us(k-1,1)*cos(Ys(k-1,3));
    Ys(k,2) = Ys(k-1,2) + Ts*Us(k-1,1)*sin(Ys(k-1,3));
    Ys(k,3) = Ys(k-1,3) + Ts*Us(k-1,2); 

    position_k = Ys(k,:);

    if position_k(3) > 2*pi % impedir que o theta ultrapasse 2pi ou -2pi
        Ys(k,3) = mod(Ys(k,3), 2*pi);
    elseif position_k(3) < -2*pi
        Ys(k,3) = mod(Ys(k,3), -2*pi);
    end

    % calculo f
    eta_x_k = Ys(k,1) - (Gxv(1,1)*dUs(k-1,1) + resposta_livre(1,1));
    eta_y_k = Ys(k,2) - (Gyv(1,1)*dUs(k-1,1) + resposta_livre(1,2));
    eta_theta_k = Ys(k,3) - (Gthetaw(1,1)*dUs(k-1,2) + resposta_livre(1,3));

    resposta_livre = modelo_com_perturbacao(Ys(k,1), Ys(k,2), Ys(k,3), Us(k-1,1), Us(k-1,2), eta_x_k, eta_y_k, eta_theta_k, N2, Ts);
    f_total = [resposta_livre(N1:N2,1); resposta_livre(N1:N2,2); resposta_livre(N1:N2,3)];

    % calculo das Gs
    [Gxv, Gxw, Gyv, Gyw] = get_Gs_f(Ys(k,1), Ys(k,2), Ys(k,3), Us(k-1,1), Us(k-1,2), e, N2, Nuv, Nuw, Ts);
  
    G_total = [  Gxv(N1:N2,:)      Gxw(N1:N2,:);
                 Gyv(N1:N2,:)      Gyw(N1:N2,:);
               Gthetav(N1:N2,:)  Gthetaw(N1:N2,:)];  

    Hqp = 2*(G_total'*Qe*G_total+Qu);
    Hqp = (Hqp+Hqp')/2;
    fqp = -2*G_total'*Qe*(Ref-f_total);

    % restrições:
    A_restr_saida = [ G_total; -G_total];

    b_restr_saida = b_aux_restr_saida + ...
        [-resposta_livre(N1:N2,1);-resposta_livre(N1:N2,2);-resposta_livre(N1:N2,3);
          resposta_livre(N1:N2,1); resposta_livre(N1:N2,2); resposta_livre(N1:N2,3)];

    b_restr_controle = b_aux_restr_controle + ...
        [-Us(k-1,1)*ones(Nuv,1); Us(k-1,1)*ones(Nuv,1);
         -Us(k-1,2)*ones(Nuw,1); Us(k-1,2)*ones(Nuw,1)];

    A = [A_restr_saida;
         A_restr_controle;
         A_restr_increm];

    b = [b_restr_saida;
         b_restr_controle;
         b_restr_increm];

    % variável de folga:

    Hw(1:end-n_epsilons,1:end-n_epsilons) = Hqp;
    fw(1:end-n_epsilons) = fqp;
    Aw(:,1:end-n_epsilons) = A;

    % solver
    Is(k,:) = quadprog(Hw,fw,Aw,b);

    % atualizando o controle
    dUs(k,1) = Is(k,1);
    dUs(k,2) = Is(k,Nuv+1);
    Us(k,1) = Us(k-1)+dUs(k,1);
    Us(k,2) = Us(k-1)+dUs(k,2);
    % [ eps_x_min, eps_y_min, eps_theta_min, eps_x_max, eps_y_max, eps_theta_max ]
    Epsilons(k,:) = Is(k,end-n_epsilons+1:end);

    k=k+1;
end
toc


%% ------------------------------ PLOTS -----------------------------------


fim_sim = k-1;

x=0:Ts:(fim_sim-ini_sim)*Ts;

figure
subplot(1,4,1)
scatter(Ys(ini_sim:fim_sim,1),Ys(ini_sim:fim_sim,2),'SizeData',14);
hold on
plot(Ys(ini_sim:fim_sim,1),Ys(ini_sim:fim_sim,2),'r');
hold on
ref_scatter = scatter(Referencia_x,Referencia_y,'g');
legend(ref_scatter,sprintf('ref x: %.2f ref y: %.2f ref theta: %.2f°',Referencia_x,Referencia_y,Referencia_theta*180/pi));
title(sprintf('Posição \nTs = %.2f s, N = %.2f s, Nuv = %.2f s, Nuw = %.2f s', Ts, N*Ts, Nuv*Ts, Nuw*Ts))
xlabel('X')
ylabel('Y')
grid on

subplot(1,4,2)
plot(x,Ys(ini_sim:fim_sim,1))
title('X')
xlabel('time [s]')
ylabel('X')
grid on

subplot(1,4,3)
plot(x,Ys(ini_sim:fim_sim,2))
title('Y')
xlabel('time [s]')
ylabel('Y')
grid on

subplot(1,4,4)
plot(x,Ys(ini_sim:fim_sim,3)*180/pi)
title('Theta')
xlabel('time [s]')
ylabel('Theta [°]')
grid on

figure
subplot(1,2,1)
plot(x,Us(ini_sim:fim_sim,1))
hold on
plot(x,dUs(ini_sim:fim_sim,1))
legend('v','dv')
title('vel. linear controle')
xlabel('time [s]')
ylabel('vel. linear')
grid on

subplot(1,2,2)
plot(x,dUs(ini_sim:fim_sim,2))
hold on
plot(x,Us(ini_sim:fim_sim,2))
legend('w','dw')
title('vel. angular controle')
xlabel('time [s]')
ylabel('vel. angular')
grid on

figure
subplot(1,4,1)
plot(x,Epsilons(ini_sim:fim_sim,:))
legend(sprintf('x max = %.2f',x_max),sprintf('y max = %.2f',y_max),sprintf('theta max = %.2f °',theta_max*180/pi), ...
       sprintf('x min = %.2f',x_min),sprintf('y min = %.2f',y_min),sprintf('theta min = %.2f °',theta_min*180/pi));
title('Variáveis de folga')
xlabel('time [s]')
ylabel('var. folga')

subplot(1,4,2)
plot(x,Ys(ini_sim:fim_sim,1))
title('X')
xlabel('time [s]')
ylabel('X')
grid on

subplot(1,4,3)
plot(x,Ys(ini_sim:fim_sim,2))
title('Y')
xlabel('time [s]')
ylabel('Y')
grid on

subplot(1,4,4)
plot(x,Ys(ini_sim:fim_sim,3)*180/pi)
title('Theta')
xlabel('time [s]')
ylabel('Theta [°]')
grid on

% -------------------------------------------------------------------------
% ANIMATION

dvs = [    Is(:,1:Nuv)     zeros(size(Is,1),N2-Nuv)]; % ajustando para o tamanho de N2
dws = [Is(:,Nuv+1:Nuv+Nuw) zeros(size(Is,1),N2-Nuw)];

prediction = modelo(pos_inicial(1), pos_inicial(2), pos_inicial(3), ...
    0, 0, dvs(1,:)', dws(1,:)', N2, Ts);

figure;
h1 = plot(prediction(1,1), prediction(1,2), 'bo','MarkerSize', 10);
hold on
h2 = plot(prediction(2:end,1), prediction(2:end,2), 'rX');
xlim([x_min-0.5 x_max+0.5]);
ylim([y_min-0.5 y_max+0.5]);
hold on;
plot(Ys(ini_sim:fim_sim,1),Ys(ini_sim:fim_sim,2), 'k--');
hold on
scatter(Referencia_x,Referencia_y,'g')
title(sprintf('PNMPC \nTs = %.2f s, N = %.2f s, Nuv = %.2f s, Nuw = %.2f s', Ts, N*Ts, Nuv*Ts, Nuw*Ts))
legend('carrinho','predições','trajeto final',sprintf('ref x: %.2f ref y: %.2f ref theta: %.2f°',Referencia_x,Referencia_y,Referencia_theta*180/pi))
xlabel('X')
ylabel('Y')

for k = ini_sim:fim_sim
    prediction = modelo(Ys(k,1), Ys(k,2), Ys(k,3), ...
    Us(k,1), Us(k,2), dvs(k,:)', dws(k,:)', N2, Ts);
    
    set(h1, 'XData', prediction(1,1), 'YData', prediction(1,2));
    set(h2, 'XData', prediction(2:end,1), 'YData', prediction(2:end,2));
    pause(Ts*2);
end


%% ----------------------------- MODELO -----------------------------------


% dv e dw tem tamanho N2x1
function prediction = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, dw, N2, Ts)
    x(1) = x_k;
    y(1) = y_k;
    theta(1) = theta_k;

    v = v_kmin1*ones(N2,1) + tril(ones(N2))*dv;
    w = w_kmin1*ones(N2,1) + tril(ones(N2))*dw;

    for k=1:N2
        x(k+1) = x(k) + Ts*v(k)*cos(theta(k));
        y(k+1) = y(k) + Ts*v(k)*sin(theta(k));
        theta(k+1) = theta(k) + Ts*w(k);
    end

    prediction(:,1) = x(2:end)';
    prediction(:,2) = y(2:end)';
    prediction(:,3) = theta(2:end)';
end

% modelo que inclue perturbações não medíveis
% dv e dw tem tamanho N2x1
function prediction = modelo_com_perturbacao(x_k, y_k, theta_k, v_kmin1, w_kmin1, eta_x_k, eta_y_k, eta_theta_k, N2, Ts)
    x(1) = x_k;
    y(1) = y_k;
    theta(1) = theta_k;

    for k=1:N2
        x(k+1) = x(k) + Ts*v_kmin1*cos(theta(k)) + eta_x_k;
        y(k+1) = y(k) + Ts*v_kmin1*sin(theta(k)) + eta_y_k;
        theta(k+1) = theta(k) + Ts*w_kmin1 + eta_theta_k;
    end

    prediction(:,1) = x(2:end)';
    prediction(:,2) = y(2:end)';
    prediction(:,3) = theta(2:end)';
end

%% ------------------------ FUNCOES AUXILIARES ----------------------------


% devido a linearidade na EDO de theta, as matrizes Gthetav e Gthetaw sã
% o triviais:
%   - Gthetav tem apenas termos nulos.
%   - Gthetaw -> gthetaw = [1*T; 2*T; ...; N2*T].
% Desse modo, as linhas que calculam essas Gs foram comentadas para otim
% izar a função.
% A função a seguir se baseia no algoritmos do PNMPC para obtenção da re
% sposta livre e matriz G
function [Gxv, Gxw, Gyv, Gyw, Gthetav, Gthetaw] = get_Gs_f(x_k, y_k, theta_k, v_kmin1, w_kmin1, e, N2, Nuv, Nuw, Ts)

% resposta livre
    dv = zeros(N2,1);
    dw = zeros(N2,1);
    resposta_livre = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, dw, N2, Ts);


% Gs relacionados a variavel de controle v
    Gxv=[];
    Gyv=[];
    Gthetav=[];
    for col=1:Nuv
        dv(col) = e;
        if col > 1
            dv(col-1) = 0;
        end
        
        resposta_incremento = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, zeros(N2,1), N2, Ts);

        Gxv = [Gxv, (resposta_incremento(:,1)-resposta_livre(:,1))/e];
        Gyv = [Gyv, (resposta_incremento(:,2)-resposta_livre(:,2))/e];
        % Gthetav = [Gthetav, (resposta_incremento(:,3)-resposta_livre(:,3))/e];
    end


% Gs relacionados a variavel de controle w
    Gxw=[];
    Gyw=[];
    Gthetaw=[];
    for col=1:Nuw
        dw(col) = e;
        if col > 1
            dw(col-1) = 0;
        end
        
        resposta_incremento = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, zeros(N2,1), dw, N2, Ts);

        Gxw = [Gxw, (resposta_incremento(:,1)-resposta_livre(:,1))/e];
        Gyw = [Gyw, (resposta_incremento(:,2)-resposta_livre(:,2))/e];
        % Gthetaw = [Gthetaw, (resposta_incremento(:,3)-resposta_livre(:,3))/e];
    end
end