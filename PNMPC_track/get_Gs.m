% devido a linearidade na EDO de theta, as matrizes Gthetav e Gthetaw sã
% o triviais:
%   - Gthetav tem apenas termos nulos.
%   - Gthetaw -> gthetaw = [1*T; 2*T; ...; N2*T].
% Desse modo, as linhas que calculam essas Gs foram comentadas para otim
% izar a função.
% A função a seguir se baseia no algoritmos do PNMPC para obtenção da re
% sposta livre e matriz G
function [Gxv, Gxw, Gyv, Gyw, Gthetav, Gthetaw] = get_Gs(x_k, y_k, theta_k, v_kmin1, w_kmin1, e, N2, Nuv, Nuw, Ts)

% resposta livre
    dv = zeros(N2,1);
    dw = zeros(N2,1);
    resposta_livre = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, dw, 0, 0, 0, N2, Ts);


% Gs relacionados a variavel de controle v
    Gxv=[];
    Gyv=[];
    Gthetav=[];
    for col=1:Nuv
        dv(col) = e;
        if col > 1
            dv(col-1) = 0;
        end
        
        resposta_incremento = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, dv, zeros(N2,1), 0, 0, 0, N2, Ts);

        Gxv = [Gxv, (resposta_incremento(:,1)-resposta_livre(:,1))/e];
        Gyv = [Gyv, (resposta_incremento(:,2)-resposta_livre(:,2))/e];
        Gthetav = [Gthetav, (resposta_incremento(:,3)-resposta_livre(:,3))/e];
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
        
        resposta_incremento = modelo(x_k, y_k, theta_k, v_kmin1, w_kmin1, zeros(N2,1), dw, 0, 0, 0, N2, Ts);

        Gxw = [Gxw, (resposta_incremento(:,1)-resposta_livre(:,1))/e];
        Gyw = [Gyw, (resposta_incremento(:,2)-resposta_livre(:,2))/e];
        Gthetaw = [Gthetaw, (resposta_incremento(:,3)-resposta_livre(:,3))/e];
    end
end