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

run('setup.m');

%% --------------------------- MONTAGEM --------------------------------


N = N2-N1+1;

Ref_x = Referencia_x(1)*ones(N,1);
Ref_y = Referencia_y(1)*ones(N,1);
Ref_theta = Referencia_theta(1)*ones(N,1);
Ref = [Ref_x; Ref_y; Ref_theta]; % transformando em vetor para fazer operações como -2*G_total'*Qe*(Ref-f_total)

dx = Referencia_x(1)-pos_inicial(1);
dy = Referencia_y(1)-pos_inicial(2);
dtheta = Referencia_theta(1)-pos_inicial(3);
if dx == 0
    dx = 1;
end
if dy == 0
    dy = 1;
end
if dtheta == 0
    dtheta = 1;
end

Qe_x = delta_x*eye(N)/(dx^2*N); 
Qe_y = delta_y*eye(N)/(dy^2*N);
Qe_theta = delta_theta*eye(N)/(dtheta^2*N);
Qe = [   Qe_x   zeros(N)  zeros(N);
       zeros(N)   Qe_y    zeros(N);
       zeros(N) zeros(N) Qe_theta]; % transformando em matriz para fazer operações como 2*(G_total'*Qe*G_total+Qu)


Qu_v = lambda_v*eye(Nuv)/(v_max^2*Nuv);
Qu_w = lambda_w*eye(Nuw)/(w_max^2*Nuw);
Qu = [     Qu_v      zeros(Nuv,Nuw);
      zeros(Nuw,Nuv)    Qu_w]; % transformando em matriz para fazer operações como 2*(G_total'*Qe*G_total+Qu)

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
aux_Aw = zeros(n_rows_Aw,Nuv+Nuw+n_epsilons);

TR = zeros(N * n_epsilons, n_epsilons); % top right (G;-G é o primeiro bloco da matriz A)
for i = 1:n_epsilons
    TR((i-1)*N+1:i*N, i) = -1; % matrix com os epsilons para somar com G
end

aux_Aw(1:size(TR,1) , end-n_epsilons+1:end) = TR;

% matriz Hw a partir de Hqp
BR = 2*diag([psi_x psi_y psi_theta psi_x psi_y psi_theta]); % bottom right

aux_Hw = zeros(Nuw+Nuv+n_epsilons);
aux_Hw(end-n_epsilons+1:end,end-n_epsilons+1:end) = BR;

% vetor fw a partir de fqp
aux_fw = zeros(Nuw+Nuv+n_epsilons,1);


%% --------------------------- SIMULAÇÃO ----------------------------------

% inicializações
resposta_livre = zeros(N2,3);
Gxv=[0];Gyv=[0];Gthetaw=[0];

ref_changes=[1]; % para plotar
ref_idx = 1; % para trocar de referencia

% -------------------------------------------------------------------------

ini_sim = N1+1;
fim_sim = ini_sim+simTim/Ts;

dUs = zeros(fim_sim,2);
Us  = zeros(fim_sim,2);
Ys  = zeros(fim_sim,3); Ys(1:ini_sim,:) = repmat(pos_inicial,ini_sim,1);
Is  = zeros(fim_sim,Nuv+Nuw+n_epsilons);
Epsilons = zeros(fim_sim,n_epsilons);

k=ini_sim;
tic
while (k<=fim_sim && ref_idx ~= length(Track_x))
    % calculando o processo
    Ys(k,1) = Ys(k-1,1) + Ts*Us(k-1,1)*cos(Ys(k-1,3));
    Ys(k,2) = Ys(k-1,2) + Ts*Us(k-1,2)*sin(Ys(k-1,3));
    Ys(k,3) = Ys(k-1,3) + Ts*Us(k-1,2); 

    % mudança de referencia
    if norm([Referencia_x(ref_idx) Referencia_y(ref_idx)]-Ys(k,1:2)) < tol
        ref_changes=[ref_changes;k]; % para plotar
        ref_idx=ref_idx+1;
        [Ref,Qe] = gerenciador_referencia(N,ref_idx,delta_x,delta_y,delta_theta,Referencia_x,Referencia_y,Referencia_theta);
    end

    % calculo f
    eta_x_k = 0;%Ys(k,1) - (Gxv(1,1)*dUs(k-1,1) + resposta_livre(1,1));
    eta_y_k = 0;%Ys(k,2) - (Gyv(1,1)*dUs(k-1,1) + resposta_livre(1,2));
    eta_theta_k = 0;%Ys(k,3) - (Gthetaw(1,1)*dUs(k-1,2) + resposta_livre(1,3));

    resposta_livre = modelo(Ys(k,1), Ys(k,2), Ys(k,3), Us(k-1,1), Us(k-1,2), zeros(N2,1), zeros(N2,1), eta_x_k, eta_y_k, eta_theta_k, N2, Ts);
    f_total = [resposta_livre(N1:N2,1); resposta_livre(N1:N2,2); resposta_livre(N1:N2,3)];

    % calculo das Gs
    [Gxv, Gxw, Gyv, Gyw, Gthetav, Gthetaw] = get_Gs(Ys(k,1), Ys(k,2), Ys(k,3), Us(k-1,1), Us(k-1,2), e, N2, Nuv, Nuw, Ts);
  
    G_total = [  Gxv(N1:N2,:)      Gxw(N1:N2,:);
                 Gyv(N1:N2,:)      Gyw(N1:N2,:);
               Gthetav(N1:N2,:)  Gthetaw(N1:N2,:)];

    Hqp = 2*(G_total'*Qe*G_total+Qu); Hqp = (Hqp+Hqp')/2;
    fqp = -2*G_total'*Qe*(Ref-f_total);

    % restrições:
    [Hw,fw,Aw,b] = gerenciador_restricao(G_total,resposta_livre,N1,N2,Nuv,Nuw,n_epsilons,Us(k-1,1),Us(k-1,2),aux_Hw,aux_fw,aux_Aw,A_restr_controle,A_restr_increm,b_aux_restr_saida,b_aux_restr_controle,b_restr_increm,Hqp,fqp);

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

fim_sim = k-1;

save("plot.mat",'Ts','fim_sim','ini_sim','N1','N2','N','Nuv','Nuw','pos_inicial','Referencia_x','Referencia_y','Referencia_theta', ...
     'x_max','y_max','theta_max','x_min','y_min','theta_min','ref_changes','Ys','Us','dUs','Is','Epsilons');