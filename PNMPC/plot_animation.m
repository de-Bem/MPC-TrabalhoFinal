close all

%% ------------------------------ PLOT ------------------------------------

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