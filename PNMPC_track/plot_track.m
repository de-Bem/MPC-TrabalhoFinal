clc
close all
clear all

speed=1; % velocidade da animação

%% ------------------------------ PLOT ------------------------------------

load('plot.mat')
load('track.mat')

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

%% --------------------------- ANIMATION ----------------------------------

dvs = [    Is(:,1:Nuv)     zeros(size(Is,1),N2-Nuv)]; % ajustando para o tamanho de N2
dws = [Is(:,Nuv+1:Nuv+Nuw) zeros(size(Is,1),N2-Nuw)];

prediction = modelo(pos_inicial(1), pos_inicial(2), pos_inicial(3), ...
    0, 0, dvs(1,:)', dws(1,:)', 0, 0, 0, N2, Ts);

figure;
h1 = plot(prediction(1,1), prediction(1,2), 'bo','MarkerSize', 10);
hold on
h2 = plot(prediction(2:end,1), prediction(2:end,2), 'rX');
xlim([x_min-0.5 x_max+0.5]);
ylim([y_min-0.5 y_max+0.5]);
hold on;
plot(Ys(ini_sim:fim_sim,1),Ys(ini_sim:fim_sim,2), 'k--');
hold on
scatter(Track_x,Track_y,'green')
title(sprintf('PNMPC track \nTs = %.2f s, N = %.2f s, Nuv = %.2f s, Nuw = %.2f s', Ts, N*Ts, Nuv*Ts, Nuw*Ts))
xlabel('X')
ylabel('Y')

i=1;
plot_ref = [Referencia_x(1) Referencia_y(1) Referencia_theta(1)];
for k = ini_sim:fim_sim
    prediction = modelo(Ys(k,1), Ys(k,2), Ys(k,3), ...
    Us(k-1,1), Us(k-1,2), dvs(k,:)', dws(k,:)', 0, 0, 0, N2, Ts);

    set(h1, 'XData', prediction(1,1), 'YData', prediction(1,2));
    set(h2, 'XData', prediction(2:end,1), 'YData', prediction(2:end,2));
    legend('carrinho', 'predições', 'trajeto final', sprintf('ref x: %.2f ref y: %.2f ref theta: %.2f°', plot_ref(1), plot_ref(2), plot_ref(3)));
    if ismember(k,ref_changes)
        i=i+1;
        plot_ref = [Referencia_x(i) Referencia_y(i) Referencia_theta(i)*180/pi];
    end

    pause(Ts/speed);
end